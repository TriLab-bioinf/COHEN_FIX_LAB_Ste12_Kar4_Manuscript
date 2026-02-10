#!/bin/bash

bwa=~/builds/bwa/bwa
sam=~/builds/samtools-1.19.2/samtools
pic=~/bin/picard.jar

ref=/Users/rogersjv/builds/bwa/genomes/sacCer3.fa

prefix=36147_Ste12_i5006_MY16679_Ste12-TAP_YPDpH3.5_Alpha-factor_XO
r1=${prefix}_R1.fastq.gz
r2=${prefix}_R2.fastq.gz


$bwa mem -t 4 -v 1 $ref $r1 $r2 | $sam sort -@4 -T "${TMPDIR:-.}" -O bam -o ${prefix}_raw.bam

java -jar $pic MarkDuplicates  INPUT=${prefix}_raw.bam OUTPUT=${prefix}_input.bam  METRICS_FILE=${prefix}.metrics.tab  REMOVE_DUPLICATES='false' \
	ASSUME_SORTED='true'  DUPLICATE_SCORING_STRATEGY='SUM_OF_BASE_QUALITIES'  \
        READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' OPTICAL_DUPLICATE_PIXEL_DISTANCE='100' \
        VALIDATION_STRINGENCY='LENIENT' TAGGING_POLICY=All QUIET=true VERBOSITY=ERROR

$sam view -o ${prefix}_FilteredBAMnew.bam -h -b -f 0x1 -F 0x404 ${prefix}_input.bam 2>&1
