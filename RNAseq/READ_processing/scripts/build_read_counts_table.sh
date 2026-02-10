#!/bin/bash

# Exit script on error
set -o errexit

source config.txt

# Create read_counts folder if it is not already there
mkdir -p ${READ_COUNTS}

# Store a list of all deduplicated bam files into the variable BAMS
BAMS=$(ls ${DEDUP_OUT}/*.sorted.dedup.bam)

module load subread

featureCounts -t exon -g gene_id -O \
    -s 2 -J -R BAM -p --ignoreDup -M \
    --fraction -G ${REF_FASTA} \
    -T 16 \
    -p --countReadPairs \
    -a ${REF_ANNOT_GTF} \
    -o ${READ_COUNTS}/read_counts_table ${BAMS}

echo ""
echo "Finished building table of read counts."
echo ""
