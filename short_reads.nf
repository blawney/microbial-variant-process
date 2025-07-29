nextflow.enable.dsl=2
nextflow.preview.output=true

include { shift_reference; 
          create_faidx; 
          create_ref_dict; 
          liftover_and_combine_vcf;
          merge_mutect_stats;
          filter_calls;
          shift_back_bam } from './modules/utils.nf'
include { align; 
          align as shifted_align } from './subworkflows/short_read_align.nf'
include { mutect2; 
          mutect2 as shifted_mutect2} from './subworkflows/mutect.nf'


workflow {

    main:

        sample_metadata_ch = Channel.fromPath(params.samplesheet)
                                 .splitCsv(header: true)

        // get the unique reference fastas. Note that we keep the "basename"
        // of the fasta file to link everything after
        unique_refs = sample_metadata_ch.map{
            row ->
                [row["ref_fasta"].split("/")[-1], row["ref_fasta"]]
        }.unique()

        // for each of the unique references, create the dict and fai files.
        // then join them into a single channel
        create_ref_dict_output = unique_refs | create_ref_dict
        create_faidx_output = unique_refs | create_faidx
        ref_materials = create_ref_dict_output.combine(create_faidx_output, by: 0)

        // Create the shifted reference to account for the circularized genome
        shifted_ref_output = unique_refs.combine(ref_materials, by: 0) | shift_reference

        // finally merge the "regular" and "shifted" reference files, shiftback, etc.
        // into a single channel keyed by the name of the reference
        all_reference_materials = ref_materials.combine(shifted_ref_output, by: 0).map{
            item ->
                x = [
                    ref_dict: item[1],
                    faidx: item[2],
                    shifted_fasta: item[3],
                    shifted_faidx: item[4],
                    shifted_dict: item[5],
                    shiftback_chain: item[6],
                    ref_intervals: item[7],
                    shifted_intervals: item[8],
                ]
                [item[0], x]
        }

        // merge the sample metadata with the corresponding reference files
        sample_metadata_ch = sample_metadata_ch.map{
            row ->
                [row["ref_fasta"].split("/")[-1], row]
        }.combine(all_reference_materials, by:0).map{
            item ->
                item[1] + item[2]
        }

        // Next we prepare to align against the "regular" and "shifted"
        // references. Since we use the same workflow for both, we have
        // to re-map the input channels to each so that the regular reference
        // files are fed to the standard alignment and the shifted reference
        // files are fed to the shifted alignment
        standard_align_input = sample_metadata_ch.map{
            item ->
                meta = [
                    readgroup_name: item.readgroup_name,
                    library_name: item.library_name,
                    platform_unit: item.platform_unit,
                    platform_name: item.platform_name,
                    sequencing_center: item.sequencing_center,
                    run_date: item.run_date
                ]
                [
                    item.sample_name,
                    item.fastq_1,
                    item.fastq_2,
                    item.ref_fasta,
                    item.ref_dict,
                    item.faidx,
                    meta
                ]
        }

        shifted_align_input = sample_metadata_ch.map{
            item ->
                meta = [
                    readgroup_name: item.readgroup_name,
                    library_name: item.library_name,
                    platform_unit: item.platform_unit,
                    platform_name: item.platform_name,
                    sequencing_center: item.sequencing_center,
                    run_date: item.run_date
                ]
                [
                    item.sample_name,
                    item.fastq_1,
                    item.fastq_2,
                    item.shifted_fasta,  // <- note the shifted reference
                    item.shifted_dict,   // <- note the shifted reference
                    item.shifted_faidx,  // <- note the shifted reference
                    meta
                ]
        }

        // Perform the alignments
        standard_align_output = standard_align_input | align
        shifted_align_output = shifted_align_input | shifted_align

        // mutect needs the BAMS plus info about the reference.
        // Similar to the alignments, we prepare similar channels
        // for passing to the regular and shifted mutect workflows
        mutect2_input = sample_metadata_ch.map{
            item ->
                [item['sample_name'], item]
        }.combine(standard_align_output, by: 0).map{
            item ->
                [
                    item[1]['sample_name'],
                    item[2],
                    item[3],
                    item[1]['ref_fasta'],
                    item[1]['faidx'],
                    item[1]['ref_dict'],
                    item[1]['ref_intervals']
                ]
        }

        shifted_mutect2_input = sample_metadata_ch.map{
            item ->
                [item['sample_name'], item]
        }.combine(shifted_align_output, by: 0).map{
            item ->
                [
                    item[1]['sample_name'],
                    item[2],
                    item[3],
                    item[1]['shifted_fasta'],      // <- note the shifted reference
                    item[1]['shifted_faidx'],      // <- note the shifted reference
                    item[1]['shifted_dict'],       // <- note the shifted reference
                    item[1]['shifted_intervals']   // <- note the shifted reference
                ]
        }

        mutect2_output = mutect2(mutect2_input, params.num_dangling_bases, params.make_output_bam)
        shifted_mutect2_output = shifted_mutect2(shifted_mutect2_input, params.num_dangling_bases, params.make_output_bam)

        // merge and reformat so it's not just a long list of outputs:
        merged_mutect = mutect2_output.combine(shifted_mutect2_output, by: 0).map {
            item ->
                x = item[1..5]
                y = item[6..10]
                [
                    item[0], 
                    [
                        orig_mutect: x, 
                        shifted_mutect: y
                    ]
                ]
        }

        // merge the sample metadata with the corresponding mutect results
        sample_metadata_ch = sample_metadata_ch.map {
                item ->
                    [item['sample_name'], item]
            }.combine(merged_mutect, by: 0).map{
                item -> item[1] + item[2]
            }

        liftover_and_combine_vcf_output = sample_metadata_ch.map{
            item ->
            [
                item['sample_name'],
                item['orig_mutect'][0],
                item['shifted_mutect'][0],
                item['ref_fasta'],
                item['faidx'],
                item['ref_dict'],
                item['shiftback_chain']
            ]
        } | liftover_and_combine_vcf

        merge_mutect_stats_output = sample_metadata_ch.map{
            item ->
            [
                item['sample_name'],
                item['orig_mutect'][2],
                item['shifted_mutect'][2]
            ]
        } | merge_mutect_stats

        filter_input = sample_metadata_ch.map{
            item -> [item['sample_name'], item]
        }.combine(liftover_and_combine_vcf_output, by:0)\
         .combine(merge_mutect_stats_output, by:0).map{
            item -> [
                item[0], 
                item[1]["ref_fasta"],
                item[1]["faidx"],
                item[1]["ref_dict"],
                item[3],
                item[4],
                item[5]
            ]
        }
        filter_output = filter_input | filter_calls

        // extract the bam's/bai's to publish:
        bams = standard_align_output.map{it[1]}
        bam_idxs = standard_align_output.map{it[2]}

    publish:
        final_vcf = filter_output.filtered_vcf
        final_vcf_idx = filter_output.filtered_vcf_idx
        liftover = liftover_and_combine_vcf_output
        bams = bams
        bam_idxs = bam_idxs

}

output {
    final_vcf {
        path 'filtered_vcf_output'
    }
    final_vcf_idx {
        path 'filtered_vcf_output'
    }
    liftover {
        path 'liftover_and_combine'
    }
    bams {
        path 'bam_files'
    }
    bam_idxs {
        path 'bam_files'
    }
}