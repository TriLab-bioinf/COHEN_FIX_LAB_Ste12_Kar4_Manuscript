#!/bin/bash
#SBATCH --cpus-per-task=8

# Exit script on error
set -o errexit

source config.txt

if [[ -z $1 ]]; then
  echo "Error: No sample name provided. Please provide a sample name as an argument when running trim_reads.sh (e.g. 'trim_reads.sh my_sample1_prefix')."
  exit 1
fi

# Set SAMPLE_PREFIX with sample name entered after command when running sbatch (e.g. "sbatch trim_reads.sh my_sample1_prefix" )
SAMPLE_PREFIX=$1

# Load required modules
module load trimgalore

# get input files
FASTQ_1=${SAMPLE_PREFIX}.R1.fastq.gz
FASTQ_2=${SAMPLE_PREFIX}.R2.fastq.gz
LOGFILE=${SAMPLE_PREFIX}_trim.log

mkdir -p ${TRIMMED_READS}

trim_galore --cores 8 --paired --basename ${SAMPLE_PREFIX} -o ${TRIMMED_READS} ${RAW_READS}/${FASTQ_1} ${RAW_READS}/${FASTQ_2} > ${TRIMMED_READS}/${LOGFILE} 2>&1

echo ""
echo "Finished trimming..."
echo ""
