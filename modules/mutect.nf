nextflow.enable.dsl=2

process run_mutect2 {

    cpus 2
    memory '8 GB'
    label 'gatk'

    input:
        val sample_name
        path input_bam
        path input_bai
        path ref_fasta
        path ref_faidx
        path ref_dict 
        path intervals
        val num_dangling_bases
        val make_output_bam

    output:
        path "${output_vcf}", emit: raw_vcf
        path "${output_vcf}.idx", emit: raw_vcf_idx
        path "${output_vcf}.stats", emit: vcf_stats
        path "${sample_name}.mutect2.bam", emit: bam_out
        path "${sample_name}.mutect2.bai", emit: bai_out

    script:
        output_vcf = "raw.vcf"
        bam_out_cmd = make_output_bam == true ? "--bam-output ${sample_name}.mutect2.bam" : ""
        """
        # File needs to exist, even if empty
        touch output.bam

        /gatk/gatk --java-options "-Xmx2000m" Mutect2 \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${output_vcf} \
            -L ${intervals} \
            ${bam_out_cmd} \
            --annotation StrandBiasBySample \
            --num-matching-bases-in-dangling-end-to-recover ${num_dangling_bases} \
            --max-reads-per-alignment-start 75 
        """
}

workflow mutect2 {

    take:
      sample_name
      input_bam
      input_bai
      ref_fasta
      ref_faidx
      ref_dict 
      intervals
      num_dangling_bases
      make_output_bam

    main:
        mutect_output = run_mutect2(
                sample_name,
                input_bam,
                input_bai,
                ref_fasta,
                ref_faidx,                                                                      
                ref_dict,
                intervals,
                num_dangling_bases,
                make_output_bam)

    emit:
        raw_vcf = mutect_output.raw_vcf
        raw_vcf_idx = mutect_output.raw_vcf_idx
        vcf_stats = mutect_output.vcf_stats
        output_bam = mutect_output.bam_out
        output_bam_idx = mutect_output.bai_out
}


workflow {

    mutect2(params.sample_name,
            params.input_bam,
            params.input_bai,
            params.ref_fasta,
            params.ref_faidx,
            params.ref_dict,
            params.intervals,
            params.num_dangling_bases,
            params.make_output_bam
    )

}