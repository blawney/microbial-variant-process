nextflow.enable.dsl=2
nextflow.preview.output=true

include { index_reference; 
          index_reference as index_shifted_reference; 
          shift_reference; 
          create_faidx; 
          create_ref_dict; 
          paired_fastq_to_unmapped_bam;
          liftover_and_combine_vcf;
          merge_mutect_stats;
          filter_calls;
          shift_back_bam } from './modules/utils.nf'
include { align; 
          align as shifted_align } from './modules/align.nf'
include { mutect2; 
          mutect2 as shifted_mutect2} from './modules/mutect.nf'


workflow process_single_sample {

    main:

    index_reference_output = index_reference(params.ref_fasta)
    create_faidx_output = create_faidx(params.ref_fasta)
    create_ref_dict_output = create_ref_dict(params.ref_fasta)

    paired_fastq_to_unmapped_bam_output = paired_fastq_to_unmapped_bam(params.sample_name, 
                                                                    params.fastq_1, 
                                                                    params.fastq_2, 
                                                                    params.readgroup_name, 
                                                                    params.library_name,
                                                                    params.platform_unit,
                                                                    params.run_date,
                                                                    params.platform_name,
                                                                    params.sequencing_center)

    shifted_reference_output = shift_reference(params.ref_fasta, 
                                                create_faidx_output.ref_faidx, 
                                                create_ref_dict_output.ref_dict)
    index_shifted_reference_output = index_shifted_reference(shifted_reference_output.shifted_ref_fasta)

    standard_align_output = align(
        params.sample_name,
        paired_fastq_to_unmapped_bam_output.ubam,
        params.fastq_1,
        params.fastq_2,
        create_ref_dict_output.ref_dict,
        params.ref_fasta,
        create_faidx_output.ref_faidx,
        index_reference_output.ref_amb,
        index_reference_output.ref_ann,
        index_reference_output.ref_bwt,
        index_reference_output.ref_pac,
        index_reference_output.ref_sa)

    shifted_align_output = shifted_align(
        params.sample_name,
        paired_fastq_to_unmapped_bam_output.ubam,
        params.fastq_1,
        params.fastq_2,
        shifted_reference_output.shifted_ref_dict,
        shifted_reference_output.shifted_ref_fasta,
        shifted_reference_output.shifted_ref_faidx,
        index_shifted_reference_output.ref_amb,
        index_shifted_reference_output.ref_ann,
        index_shifted_reference_output.ref_bwt,
        index_shifted_reference_output.ref_pac,
        index_shifted_reference_output.ref_sa)

    mutect2_output = mutect2(params.sample_name,
        standard_align_output.bam,
        standard_align_output.bam_idx,
        params.ref_fasta,
        create_faidx_output.ref_faidx,
        create_ref_dict_output.ref_dict,
        shifted_reference_output.unshifted_intervals,
        params.num_dangling_bases,
        params.make_output_bam)

    shifted_mutect2_output = shifted_mutect2(params.sample_name,
        shifted_align_output.bam,
        shifted_align_output.bam_idx,
        shifted_reference_output.shifted_ref_fasta,
        shifted_reference_output.shifted_ref_faidx,
        shifted_reference_output.shifted_ref_dict,
        shifted_reference_output.shifted_intervals,
        params.num_dangling_bases,
        params.make_output_bam)

    if(params.make_output_bam){
        shift_back_bam_output = shift_back_bam(
            shifted_mutect2_output.output_bam,
            shifted_reference_output.shiftback_chain
        )
    }

    liftover_and_combine_vcf_output = liftover_and_combine_vcf(
        mutect2_output.raw_vcf,
        shifted_mutect2_output.raw_vcf,
        params.ref_fasta,
        create_faidx_output.ref_faidx,
        create_ref_dict_output.ref_dict,
        shifted_reference_output.shiftback_chain,
        params.sample_name
    )

    merge_mutect_stats_output = merge_mutect_stats(
        mutect2_output.vcf_stats,
        shifted_mutect2_output.vcf_stats
    )

    filter_output = filter_calls(
        params.sample_name,
        params.ref_fasta,
        create_faidx_output.ref_faidx,
        create_ref_dict_output.ref_dict,
        liftover_and_combine_vcf_output.merged_vcf,
        liftover_and_combine_vcf_output.merged_vcf_idx,
        merge_mutect_stats_output.stats
    )

    emit:
    merged_vcf = liftover_and_combine_vcf_output.merged_vcf
    merged_vcf_idx = liftover_and_combine_vcf_output.merged_vcf_idx
    final_vcf = filter_output.filtered_vcf
    final_vcf_idx = filter_output.filtered_vcf_idx
    orig_bam = standard_align_output.bam
    orig_bam_idx = standard_align_output.bam_idx
    mutect_bam = mutect2_output.output_bam
    mutect_bam_idx = mutect2_output.output_bam_idx
    
}


workflow {

    outputs = process_single_sample()

    publish:
    merged_vcf = outputs.merged_vcf
    merged_vcf_idx = outputs.merged_vcf_idx
    final_vcf = outputs.final_vcf
    final_vcf_idx = outputs.final_vcf_idx
    orig_bam = outputs.orig_bam
    orig_bam_idx = outputs.orig_bam_idx
    mutect_bam = outputs.mutect_bam
    mutect_bam_idx = outputs.mutect_bam_idx

}

output {
    merged_vcf {
        path 'unfiltered_vcf_output'
    }
    merged_vcf_idx {
        path 'unfiltered_vcf_output'
    }
    final_vcf {
        path 'filtered_vcf_output'
    }
    final_vcf_idx {
        path 'filtered_vcf_output'
    }
    orig_bam {
        path 'original_bams'
    }
    orig_bam_idx {
        path 'original_bams'
    }
    mutect_bam {
        path 'mutect2_bams'
    }
    mutect_bam_idx {
        path 'mutect2_bams'
    }
}