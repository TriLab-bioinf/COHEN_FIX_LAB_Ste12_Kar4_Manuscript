
# Scripts for processing raw paired-end RNA-seq reads

The scripts in this directory were used to build a table containing read counts per sample and per gene (samples as columns, genes as rows).

Before running the scripts, edit the `config.txt` file to set the paths for: the folder with raw FASTQ files, input and output directories for each step, the STAR index folder for the reference genome, the reference genome FASTA, and the reference annotation GTF.

Run the scripts in the following order:

```bash
# Run FastQC on raw FASTQ files
fastqc_raw.sh

# Trim raw reads with Trim Galore
trim_reads.sh

# Run FastQC on trimmed FASTQ files
fastqc_trimmed.sh

# Map reads to the reference genome with STAR
map_reads.sh

# Mark PCR duplicates using Picard
rm_duplicates.sh

# Generate a read-count table (genes x samples)
build_read_counts_table.sh
```

The scripts use `module load` to load bioinformatics tools. If `module` is not available on your system, ensure the following programs are available in one of the directories listed in your `$PATH` environment variable:

#### Required tools:

```bash
fastqc
trimgalore
STAR
picard
subread
```

Raw FASTQ file names are expected to follow this format:

```bash
${sample_id}.R1.fastq.gz
${sample_id}.R2.fastq.gz
```

`${sample_id}` is the sample name. For example, for a sample named `my_sample1`:

```bash
my_sample1.R1.fastq.gz
my_sample1.R2.fastq.gz
```
**NOTE:** These scripts provide the explicit commands used to perform the analyses described in the paper. However, they are not intended to be functional for all users without modification. Users are expected to install the required programs on their own machine (or computing cluster), then update the files to point at their own filepaths.
