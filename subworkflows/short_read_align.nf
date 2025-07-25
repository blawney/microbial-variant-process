nextflow.enable.dsl=2


process get_bwa_version {

    label 'bwa_and_picard'

    output:
        stdout

    script:
    """
    /usr/local/bin/bwa 2>&1 | grep -e '^Version' | sed 's/Version: //'
    """
}


process index_reference {

    cpus 2
    memory '4 GB'
    label 'bwa_and_picard'

    input:
        tuple val(fasta_id), path(ref_fasta)

    output:
        tuple val(fasta_id), \
            path("${ref_fasta_base}.amb"), \
            path("${ref_fasta_base}.ann"), \
            path("${ref_fasta_base}.bwt"), \
            path("${ref_fasta_base}.pac"), \
            path("${ref_fasta_base}.sa")

    script:
        ref_fasta_base = ref_fasta.name
        """
        /usr/local/bin/bwa index ${ref_fasta_base}
        """
}


process paired_fastq_to_unmapped_bam {

    cpus 2
    memory '12 GB'
    label 'gatk'

    input:
        tuple val(sample_name), \
            path(fastq_1), \
            path(fastq_2), \
            path(ref_fasta), \
            path(ref_dict), \
            path(ref_faidx), \
            val(sequencing_metadata)

    output:
        tuple val(sample_name), path("${sample_name}.unmapped.bam")

    script:
        sequencing_center = sequencing_metadata["sequencing_center"]
        readgroup_name = sequencing_metadata["readgroup_name"]
        library_name = sequencing_metadata["library_name"]
        platform_unit = sequencing_metadata["platform_unit"]
        platform_name = sequencing_metadata["platform_name"]
        run_date = sequencing_metadata["run_date"]
        """
        /gatk/gatk --java-options "-Xmx10g" \
            FastqToSam \
            --FASTQ ${fastq_1} \
            --FASTQ2 ${fastq_2} \
            --OUTPUT ${sample_name}.unmapped.bam \
            --READ_GROUP_NAME ${readgroup_name} \
            --SAMPLE_NAME ${sample_name} \
            --LIBRARY_NAME ${library_name} \
            --PLATFORM_UNIT ${platform_unit} \
            --RUN_DATE ${run_date} \
            --PLATFORM ${platform_name} \
            --SEQUENCING_CENTER ${sequencing_center} 
            """
}


process align_and_mark_duplicates {

    cpus 2
    memory '8 GB'
    label 'bwa_and_picard'

    input:
        tuple val(sample_name), \
            path(fastq_1), \
            path(fastq_2), \
            path(ref_fasta), \
            path(ref_dict), \
            path(faidx), \
            val(meta), \
            path(unmapped_bam), \
            path(ref_amb), \
            path(ref_ann), \
            path(ref_bwt), \
            path(ref_pac), \
            path(ref_sa)
        val bwa_version

    output:
        tuple val(sample_name), \
            path("${output_bam_basename}.bam"), \
            path("${output_bam_basename}.bai"), \
            path("${output_bam_basename}.bwa.stderr.log"), \
            path("${metrics_filename}")

    script:

        bwa_version_strip = bwa_version.replaceAll(/\n/, '')
        basename = "${sample_name}"
        output_bam_basename = "${basename}.realigned"
        ref_fasta_base = ref_fasta.baseName
        metrics_filename = "${basename}.metrics"
        bwa_call = "bwa mem -K 100000000 -v 3 -t 2"
        """
       /usr/local/bin/${bwa_call} \
           ${ref_fasta} \
           ${fastq_1} \
           ${fastq_2} 2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) > ${basename}.realigned_bam

        java -Xms2000m -jar /usr/local/bin/picard.jar \
            MergeBamAlignment \
            VALIDATION_STRINGENCY=SILENT \
            EXPECTED_ORIENTATIONS=FR \
            ATTRIBUTES_TO_RETAIN=X0 \
            ATTRIBUTES_TO_REMOVE=NM \
            ATTRIBUTES_TO_REMOVE=MD \
            ALIGNED_BAM=${basename}.realigned_bam \
            UNMAPPED_BAM=${unmapped_bam} \
            OUTPUT=mba.bam \
            REFERENCE_SEQUENCE=${ref_fasta} \
            PAIRED_RUN=true \
            SORT_ORDER="unsorted" \
            IS_BISULFITE_SEQUENCE=false \
            ALIGNED_READS_ONLY=false \
            CLIP_ADAPTERS=false \
            MAX_RECORDS_IN_RAM=2000000 \
            ADD_MATE_CIGAR=true \
            MAX_INSERTIONS_OR_DELETIONS=-1 \
            PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
            PROGRAM_RECORD_ID="bwamem" \
            PROGRAM_GROUP_VERSION="${bwa_version_strip}" \
            PROGRAM_GROUP_COMMAND_LINE="${bwa_call}" \
            PROGRAM_GROUP_NAME="bwamem" \
            UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
            ALIGNER_PROPER_PAIR_FLAGS=true \
            UNMAP_CONTAMINANT_READS=true \
            ADD_PG_TAG_TO_READS=false

        java -Xms2000m -jar /usr/local/bin/picard.jar \
            MarkDuplicates \
            INPUT=mba.bam \
            OUTPUT=md.bam \
            METRICS_FILE=${metrics_filename} \
            VALIDATION_STRINGENCY=SILENT \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
            ASSUME_SORT_ORDER="queryname" \
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
        sample_metadata

    main:

        bwa_version = get_bwa_version()

        // create the index for each unique fasta. Keep a reference
        // to the reference fasta so we can use that to join. This way
        // we are not creating the same index over and over for multiple
        // samples aligning to the same reference
        unique_refs = sample_metadata.map{
            row ->
                [row[3].toString().split("/")[-1], row[3]]
        }.unique()
        idx_ch = index_reference(unique_refs)

        // This channel returns a tuple of the sample name and BAM.
        ubam_output = paired_fastq_to_unmapped_bam(sample_metadata)

        // using the sample name, join the metadata to the unmapped BAM:
        sample_metadata = sample_metadata.combine(ubam_output, by: 0)

        // now need to link up the BWA index (which is obviously at
        // the reference fasta level) with the sample level info
        aligner_input_ch = sample_metadata.map{
            row ->
                [row[3].toString().split("/")[-1], row]
        }.combine(idx_ch, by:0).map{
            item ->
                // note that this is intentional. item[1]
                // is itself a list and when we add item[2..6]
                // we are appending to that. This way the sample 
                // name and fastqs aren't in a nested list
                item[1] + item[2..6]
        }

        align_out = align_and_mark_duplicates(aligner_input_ch, bwa_version)

    emit:
        alignment_outputs = align_out 
}

workflow {

    sample_metadata_ch = Channel.fromPath(params.samplesheet)
                                .splitCsv(header: true)
    align(sample_metadata_ch)
}