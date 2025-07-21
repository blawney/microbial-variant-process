nextflow.enable.dsl=2
nextflow.preview.output=true

include { process_single_sample } from './subworkflows/single_sample_process.nf'


workflow {

    main:
        sample_metadata = Channel.fromPath(params.samplesheet)
                                 .splitCsv(header: true)
                                 .map{ row ->
                                    meta = [
                                        readgroup_name: row.readgroup_name,
                                        library_name: row.library_name,
                                        platform_unit: row.platform_unit,
                                        platform_name: row.platform_name,
                                        sequencing_center: row.sequencing_center,
                                        run_date: row.run_date
                                    ]
                                    [
                                        ref_fasta: row.ref_fasta, 
                                        sample_name: row.sample_name,
                                        fastq_1: row.fastq_1,
                                        fastq_2: row.fastq_2,
                                        sequencing_metadata: meta
                                    ]
                                 }

        outputs = process_single_sample(
            sample_metadata.map{it['ref_fasta']}, 
            sample_metadata.map{it['sample_name']}, 
            sample_metadata.map{it['fastq_1']},
            sample_metadata.map{it['fastq_2']},
            sample_metadata.map{it['sequencing_metadata']},
            params.num_dangling_bases,
            params.make_output_bam
        )

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