#!/bin/bash

# Exit script on error
set -o errexit

source config.txt

if [[ -z $1 ]]; then
  echo "Error: No sample name provided. Please provide a sample name as an argument when running rm_duplicates.sh (e.g. 'rm_duplicates.sh my_sample1_prefix')."
  exit 1
fi

SAMPLE_PREFIX=$1

# Create dedup folder if it is not already there
mkdir -p ${DEDUP_OUT}

module load picard

# Create a new BAM file with duplicated reads marked (but not removed) and a metrics file with duplication statistics
# Set PICARDJAR variable in config.txt to the path of the Picard jar file (e.g. /path/to/picard.jar) before running rm_duplicates.sh if it is not already set with 'module load picard' command
# Set java memory (-Xmx) according to available resources (e.g. -Xmx32g for 32 GB of memory)
java -Xmx32g -jar ${PICARDJAR} MarkDuplicates \
    I=${STAR_OUT}/${SAMPLE_PREFIX}.Aligned.sortedByCoord.out.bam \
    O=${DEDUP_OUT}/${SAMPLE_PREFIX}.sorted.dedup.bam \
    M=${DEDUP_OUT}/${SAMPLE_PREFIX}.sorted.dedup.metrics.txt \
    READ_NAME_REGEX=null \
    REMOVE_DUPLICATES=false

echo ""
echo "Finished deduplicating reads."
echo ""