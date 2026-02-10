#!/bin/bash

# Exit script on error
set -o errexit

source config.txt

# Load fastqc tool
module load fastqc

# Create directory to store fastqc output of raw fastq files
mkdir -p -m 777 fastqc_raw

# Save path to read fastq files into the variable READS
READS=$(ls ./${TRIMMED_READS}/*.fastq.gz)

# Run fastqc on raw fastq files
fastqc --quiet -t 8 -o ${TRIMMED_READS} ${READS}
