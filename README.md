# Microbial variant calling process

This repository contains a couple Nextflow-based processes to perform alignment and variant calling using Mutect2. Based on the WDL pipelines here: [GATK-for-Microbes](https://github.com/broadinstitute/GATK-for-Microbes/tree/master/wdl/shortReads)

**All processes here assume we are working with a circularized genome**. The process performs alignment and variant-calling on both the "original" linearized reference FASTA and a shifted FASTA to account for potential alignments that span the artificial breakpoint.

There are two versions of the process for short and long reads, respectively. Note that each is able to handle situations where there are potentially multiple genomes. For instance, if you are aligning samples A,B, and C to reference genomes X, Y, and Z, it is capable of handling this.

Note that this process is quite general and there is no out-of-the-box consideration for skipping variant calling in particular regions (e.g. skipping genes associated with antibiotic resistance, etc.) 

### Short read version

The short read version assumes the following:
- Paired-end reads amenable to alignment with BWA-mem.
- Input data in FASTQ format (separate R1 and R2 files)

**To use**:

Create a metadata file in CSV-format with the following fields. The ordering does not matter, but the name of the fields in the header line *does matter*.
- `sample_name`*: Unique sample ID
- `ref_fasta`*: Path to a reference FASTA-format file
- `fastq_1`*: Path to R1
- `fastq_2`*: Path to R2
- `readgroup_name`*: Readgroup. Required by GATK tools 
- `library_name`
- `platform_name`
- `platform_unit`
- `sequencing_center`
- `run_date`

Many of the fields (those *not* marked with stars) can be left blank without causing issues with GATK.

Then run:
```
nextflow run short_reads.nf \
  -c nextflow.config \
  --samplesheet=<PATH TO METADATA CSV>
```

### Long read version

The long read version assumes the following: 
- PacBio reads which have already been run through a consensus sequence process (e.g. CCS or other)
- Input is provided in unaligned BAM format (with caveats on the BAM headers, as noted below)

**To use**:

Create a metadata file in CSV-format with the following fields. The ordering does not matter, but the name of the fields in the header line *does matter*.
- `sample_name`*: Unique sample ID
- `ref_fasta`*: Path to a reference FASTA-format file
- `unaligned_bam`*: Path to the unaligned BAM file.

Note that this assumes the unaligned BAM already has many of the required BAM headers required by GATK (e.g. `@RG`/read groups)

Then run:
```
nextflow run long_reads.nf \
  -c nextflow.config \
  --samplesheet=<PATH TO METADATA CSV>
```

### Output files

Both processes deposit files in the `results/` folder within a timestamped subdirectory, e.g. `results/YYYYMMDD-HH-mm-ss/`. Within that, you will see subdirectories for BAM files, raw VCF files, and filtered VCF files. 