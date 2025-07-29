nextflow.enable.dsl=2


process index_reference {

    cpus 2
    memory '4 GB'
    label 'longread'

    input:
        tuple val(fasta_id), path(ref_fasta)

    output:
        tuple val(fasta_id), \
            path("${ref_fasta_base}.mmi")

    script:
        ref_fasta_base = ref_fasta.name
        """
        pbmm2 index --preset CCS ${ref_fasta} ${ref_fasta_base}.mmi
        """
}


process align_and_mark_duplicates {

    cpus 2
    memory '8 GB'
    label 'longread'

    input:
        tuple val(sample_name), \
            path(unaligned_bam), \
            path(ref_fasta), \
            path(idx_mmi)

    output:
        tuple val(sample_name), \
            path("${output_bam_basename}.bam"), \
            path("${output_bam_basename}.bai"), \
            path("${metrics_filename}")

    script:

        basename = "${sample_name}"
        output_bam_basename = "${basename}.aligned_and_md"
        ref_fasta_base = ref_fasta.baseName
        metrics_filename = "${basename}.metrics"
        """
        /opt/conda/bin/pbmm2 align --sort --preset CCS \
            ${idx_mmi} ${unaligned_bam} ${basename}.init_aligned.bam

        java -Xms2000m -jar /usr/local/bin/picard.jar \
            MarkDuplicates \
            INPUT=${basename}.init_aligned.bam \
            OUTPUT=md.bam \
            METRICS_FILE=${metrics_filename} \
            VALIDATION_STRINGENCY=SILENT \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
            ASSUME_SORTED=true \
            CLEAR_DT=false \
            ADD_PG_TAG_TO_READS=false

        java -Xms2000m -jar /usr/local/bin/picard.jar \
            SortSam \
            INPUT=md.bam \
            OUTPUT=${output_bam_basename}.bam \
            SORT_ORDER="coordinate" \
            CREATE_INDEX=true \
            MAX_RECORDS_IN_RAM=300000
        """
}


workflow align {

    take:
        // note that sample_metadata looks like:
        // sample_name,
        // unaligned_bam,
        // ref_fasta
        sample_metadata

    main:

        // create the index for each unique fasta. Keep a reference
        // to the reference fasta so we can use that to join. This way
        // we are not creating the same index over and over for multiple
        // samples aligning to the same reference
        unique_refs = sample_metadata.map{
            row ->
                [row[2].toString().split("/")[-1], row[2]]
        }.unique()

        // returns a tuple of the fasta name and the MMI index
        idx_ch = index_reference(unique_refs)

        // now need to link up the pbmm2 index (which is obviously at
        // the reference fasta level) with the sample level info
        aligner_input_ch = sample_metadata.map{
            row ->
                [row[2].toString().split("/")[-1], row]
        }.combine(idx_ch, by:0).map{
            item ->
                // note that this is intentional. item[1]
                // is itself a list and when we add item[2..6]
                // we are appending to that. This way the sample 
                // name and fastqs aren't in a nested list
                item[1] + [item[2]]
        }
        aligner_input_ch.view()
        align_out = align_and_mark_duplicates(aligner_input_ch)

    emit:
        alignment_outputs = align_out 
}

workflow {

    sample_metadata_ch = Channel.fromPath(params.samplesheet)
                                .splitCsv(header: true)
    align(sample_metadata_ch)
}