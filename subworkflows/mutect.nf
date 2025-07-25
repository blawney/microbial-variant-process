nextflow.enable.dsl=2


process run_mutect2 {

    cpus 2
    memory '8 GB'
    label 'gatk'

    input:
        tuple val(sample_name), \
            path(input_bam), \
            path(input_bai), \
            path(ref_fasta), \
            path(ref_faidx), \
            path(ref_dict ), \
            path(intervals)
        val num_dangling_bases
        val make_output_bam

    output:
        tuple val(sample_name), \
            path("${output_vcf}"), \
            path("${output_vcf}.idx"), \
            path("${output_vcf}.stats"), \
            path("${sample_name}.mutect2.bam"), \
            path("${sample_name}.mutect2.bai")

    script:
        output_vcf = "${sample_name}.raw.vcf"
        bam_out_cmd = make_output_bam == true ? "--bam-output ${sample_name}.mutect2.bam" : ""
        """
        # File needs to exist, even if empty
        touch "${sample_name}.mutect2.bam"

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
        sample_metadata
        num_dangling_bases
        make_output_bam

    main:
        mutect_output = run_mutect2(
                sample_metadata,
                num_dangling_bases,
                make_output_bam)

    emit:
        mutect_output
}


workflow {

    sample_metadata_ch = Channel.fromPath(params.samplesheet)
                                .splitCsv(header: true)
                                .map {
                                    row ->
                                    [
                                        row.sample_name,
                                        row.input_bam,
                                        row.input_bai,
                                        row.ref_fasta,
                                        row.ref_faidx,
                                        row.ref_dict,
                                        row.intervals
                                    ]
                                }
    mutect2(sample_metadata_ch
            params.num_dangling_bases,
            params.make_output_bam
    )

}