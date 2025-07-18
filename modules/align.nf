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

process align_and_mark_duplicates {

    cpus 2
    memory '8 GB'
    label 'bwa_and_picard'

    input:
        path unmapped_bam
        path fastq_1
        path fastq_2
        val bwa_version
        path ref_dict
        path ref_fasta
        path ref_fasta_index
        path ref_amb
        path ref_ann
        path ref_bwt
        path ref_pac
        path ref_sa

    output:
        path "${output_bam_basename}.bam", emit: bam 
        path "${output_bam_basename}.bai", emit: bam_idx
        path "${output_bam_basename}.bwa.stderr.log", emit: bwa_stderr
        path "${metrics_filename}", emit: duplicate_metrics

    script:
        bwa_version_strip = bwa_version.replaceAll(/\n/, '')
        basename = unmapped_bam.baseName
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
        unmapped_bam
        fastq_1
        fastq_2
        ref_dict
        ref_fasta
        ref_fasta_index
        ref_amb
        ref_ann
        ref_bwt
        ref_pac
        ref_sa

    main:
        bwa_version = get_bwa_version()
        align_out = align_and_mark_duplicates(unmapped_bam,
                                    fastq_1,
                                    fastq_2,
                                    bwa_version,
                                    ref_dict,
                                    ref_fasta,
                                    ref_fasta_index,
                                    ref_amb,
                                    ref_ann,
                                    ref_bwt,
                                    ref_pac,
                                    ref_sa)

    emit:
        bam = align_out.bam 
        bam_idx = align_out.bam_idx
        align_stderr = align_out.bwa_stderr
        align_metrics = align_out.duplicate_metrics
}

workflow {
    align(params.unmapped_bam,
        params.fastq_1,
        params.fastq_2,
        params.ref_dict,
        params.ref_fasta,
        params.ref_fasta_index,
        params.ref_amb,
        params.ref_ann,
        params.ref_bwt,
        params.ref_pac,
        params.ref_sa)
}