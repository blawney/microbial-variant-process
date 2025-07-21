process shift_reference {

    cpus 2
    memory '4 GB'
    label 'gatk'

    input:
        path ref_fasta
        path ref_fasta_index
        path ref_dict

    output:
        path "${ref_fasta_base}.shifted.fa", emit: shifted_ref_fasta
        path "${ref_fasta_base}.shifted.fa.fai", emit: shifted_ref_faidx
        path "${ref_fasta_base}.shifted.dict", emit: shifted_ref_dict
        path "${ref_fasta_base}.shiftback.chain", emit: shiftback_chain
        path "${ref_fasta_base}.intervals", emit: unshifted_intervals
        path "${ref_fasta_base}.shifted.intervals", emit: shifted_intervals

    script:
        ref_fasta_base = ref_fasta.baseName
        """
        gatk --java-options "-Xmx2500m" ShiftFasta \
            -R ${ref_fasta_base}.fa \
            -O ${ref_fasta_base}.shifted.fa \
            --interval-file-name ${ref_fasta_base} \
            --shift-back-output ${ref_fasta_base}.shiftback.chain
        """
}


process index_reference {

    cpus 2
    memory '4 GB'
    label 'bwa_and_picard'

    input:
        path ref_fasta

    output:
        path "${ref_fasta_base}.amb", emit: ref_amb
        path "${ref_fasta_base}.ann", emit: ref_ann
        path "${ref_fasta_base}.bwt", emit: ref_bwt
        path "${ref_fasta_base}.pac", emit: ref_pac
        path "${ref_fasta_base}.sa", emit: ref_sa

    script:
        ref_fasta_base = ref_fasta.name
        """
        /usr/local/bin/bwa index ${ref_fasta_base}
        """

}


process create_ref_dict {

    cpus 2
    memory '4 GB'
    label 'bwa_and_picard'

    input:
        path ref_fasta

    output:
        path "${ref_fasta_prefix}.dict", emit: ref_dict

    script:
        ref_fasta_prefix = ref_fasta.baseName
        """
        java -jar /usr/local/bin/picard.jar CreateSequenceDictionary \
            R=${ref_fasta} \
            O=${ref_fasta_prefix}.dict
        """
}


process create_faidx {

    cpus 2
    memory '4 GB'
    container 'docker.io/biocontainers/samtools:v1.9-4-deb_cv1'

    input:
        path ref_fasta

    output:
        path "${ref_fasta}.fai", emit: ref_faidx

    script:
        """
        samtools faidx ${ref_fasta}
        """
}


process paired_fastq_to_unmapped_bam {

    cpus 2
    memory '12 GB'
    label 'gatk'

    input:
        val sample_name
        path fastq_1
        path fastq_2
        val sequencing_metadata

    output:
        path "${sample_name}.unmapped.bam", emit: ubam

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


process liftover_and_combine_vcf {

    cpus 2
    memory '10 GB'
    label 'bwa_and_picard'

    input:
        // these are effectively aliases, or else you get name collisions since 
        // both mutect2 processes create VCFs with the same name.
        path "unshifted_vcf"
        path "shifted_vcf"
        path ref_fasta
        path ref_faidx
        path ref_dict
        path shiftback_chain
        val sample_name

    output:
        path "${sample_name}.rejected.vcf", emit: rejected_vcf
        path "${sample_name}.merged.vcf", emit: merged_vcf
        path "${sample_name}.merged.vcf.idx", emit: merged_vcf_idx

    script:
        """
        # note that the links are necessary since picard uses the ".vcf"
        # suffix to infer file formats.
        ln -s ${shifted_vcf} shifted.vcf
        ln -s ${unshifted_vcf} unshifted.vcf

        java -jar /usr/local/bin/picard.jar LiftoverVcf \
            I=shifted.vcf \
            O=${sample_name}.shifted_back.vcf \
            R=${ref_fasta} \
            CHAIN=${shiftback_chain} \
            REJECT=${sample_name}.rejected.vcf

        java -jar /usr/local/bin/picard.jar MergeVcfs \
            I=${sample_name}.shifted_back.vcf \
            I=unshifted.vcf \
            O=${sample_name}.merged.vcf
        """     
}

process merge_mutect_stats {
    cpus 1
    memory '2 GB'
    label 'gatk'

    input:
        path "unshifted_stats"
        path "shifted_stats"

    output:
        path 'raw.combined.stats', emit: stats

    script:
        """
        /gatk/gatk MergeMutectStats \
            --stats ${shifted_stats} \
            --stats ${unshifted_stats} \
            -O raw.combined.stats
        """
}

process filter_calls {
    cpus 2
    memory '4 GB'
    label 'gatk'

    input:
        val sample_name
        path ref_fasta
        path ref_faidx
        path ref_dict
        path raw_vcf
        path raw_vcf_idx
        path raw_vcf_stats

    output:
        path "${output_vcf}", emit: filtered_vcf
        path "${output_vcf}.idx", emit: filtered_vcf_idx

    script:
        output_vcf = "${sample_name}.filtered.vcf"
        """
        /gatk/gatk --java-options "-Xmx3000m" FilterMutectCalls \
            -V ${raw_vcf} \
            -R ${ref_fasta} \
            -O ${output_vcf} \
            --stats ${raw_vcf_stats} \
            --microbial-mode 
        """
}

process shift_back_bam {
    cpus 2
    memory '8 GB'
    label 'crossmap'

    input:
        path bamfile
        path shiftback_chain

    output:
        path "bamout.sorted.bam", emit: backshifted_bam
        path "bamout.sorted.bam.bai", emit: backshifted_bam_idx

    script:
        """
        CrossMap bam ${shiftback_chain} ${bamfile} bamout
        """
}