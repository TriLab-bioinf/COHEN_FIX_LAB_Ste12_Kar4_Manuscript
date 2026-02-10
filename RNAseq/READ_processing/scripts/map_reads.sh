#!/bin/bash

# Exit script on error
set -o errexit

source config.txt

if [[ -z $1 ]]; then
  echo "Error: No sample name provided. Please provide a sample name as an argument when running map_reads.sh (e.g. 'map_reads.sh my_sample1_prefix')."
  exit 1
fi

# Check that STAR indexes exist for the reference genome
if [[ ! -d $GENOME ]]; then
  echo "Error: '${GENOME}' does not exist or is not a directory. Please set the path to the STAR genome index in config.txt before running map_reads.sh."
  exit 1
fi

# Set SAMPLE_PREFIX with sample name entered after command when running sbatch (e.g. "sbatch trim_reads.sh my_sample1_prefix" )
SAMPLE_PREFIX=$1

module load STAR

# Create mapped_reads folder if it is not already there
mkdir -p ${STAR_OUT}

STAR --runMode alignReads \
        --runThreadN 16 \
        --genomeDir ${GENOME} \
        --alignSJDBoverhangMin 1 \
        --alignSJoverhangMin 5 \
        --outFilterMismatchNmax 2 \
        --alignEndsType EndToEnd \
        --readFilesIn ${TRIMMED_READS}/${SAMPLE_PREFIX}_val_1.fq.gz ${TRIMMED_READS}/${SAMPLE_PREFIX}_val_2.fq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix ${STAR_OUT}/${SAMPLE_PREFIX}. \
        --quantMode GeneCounts \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes All

echo ""
echo "Finished mapping..."
echo ""
