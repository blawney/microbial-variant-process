process shift_reference {

    cpus 2
    memory '4 GB'
    label 'gatk'

    input:
        tuple val(ref_fasta_id), path(ref_fasta), path(ref_dict), path(ref_fasta_index)

    output:
        tuple val(ref_fasta_id), \
        path("${ref_fasta_base}.shifted.fa"), \
        path("${ref_fasta_base}.shifted.fa.fai"), \
        path("${ref_fasta_base}.shifted.dict"), \
        path("${ref_fasta_base}.shiftback.chain"), \
        path("${ref_fasta_base}.intervals"), \
        path("${ref_fasta_base}.shifted.intervals")

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


process create_ref_dict {

    cpus 2
    memory '4 GB'
    label 'bwa_and_picard'

    input:
        tuple val(ref_fasta_id), path(ref_fasta)

    output:
        tuple val(ref_fasta_id), path("${ref_fasta_prefix}.dict")

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
        tuple val(ref_fasta_id), path(ref_fasta)

    output:
        tuple val(ref_fasta_id), path("${ref_fasta}.fai")

    script:
        """
        samtools faidx ${ref_fasta}
        """
}


process liftover_and_combine_vcf {

    cpus 2
    memory '10 GB'
    label 'bwa_and_picard'

    input:
        // these are effectively aliases, or else you get name collisions since 
        // both mutect2 processes create VCFs with the same name.
        tuple val(sample_name), \
            path("unshifted_vcf"), \
            path("shifted_vcf"), \
            path(ref_fasta), \
            path(ref_faidx), \
            path(ref_dict), \
            path(shiftback_chain)

    output:
        tuple val(sample_name), \
            path("${sample_name}.rejected.vcf"), \
            path("${sample_name}.merged.vcf"), \
            path("${sample_name}.merged.vcf.idx")

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
        tuple val(sample_name), \
            path("unshifted_stats"), \
            path("shifted_stats")

    output:
        tuple val(sample_name), path("${sample_name}.combined.stats")

    script:
        """
        /gatk/gatk MergeMutectStats \
            --stats ${shifted_stats} \
            --stats ${unshifted_stats} \
            -O ${sample_name}.combined.stats
        """
}

process filter_calls {
    cpus 2
    memory '4 GB'
    label 'gatk'

    input:
        tuple val(sample_name), \
            path(ref_fasta), \
            path(ref_faidx), \
            path(ref_dict), \
            path(raw_vcf), \
            path(raw_vcf_idx), \
            path(raw_vcf_stats)

    output:
        path("${output_vcf}"), emit: filtered_vcf
        path("${output_vcf}.idx"), emit: filtered_vcf_idx

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