
#!/bin/bash

sam=~/builds/samtools-1.19.2/samtools

#prefix=36203_Ste12_i5006_MY16679_Ste12-TAP_YPDpH3.5_Alpha-factor_XO.sorted.filtered
#prefix=36205_Ste12_i5006_MY16682kar4--kanMX_Ste12-TAP_YPDpH3.5_Alpha-factor_XO_FilteredBAM

for i in *FilteredBAM.bam
do
	prefix=${i/.bam/}
	echo $prefix

	#clean bams with proper pair flag and MQ 30
	$sam view -f 2 -b -q 30 ${prefix}.bam > ${prefix}.proper.q30.bam
	$sam index ${prefix}.proper.q30.bam

	#for investigating filtered reads, keeps only improper paired reads
	#$sam view -F 2 -b ${prefix}.bam > ${prefix}.improper.bam

	#make tag count bigwigs for new bams
	~/.local/bin/bamCoverage --numberOfProcessors max/2 --bam ${prefix}.proper.q30.bam --outFileName ${prefix}.proper.q30_top.bw --outFileFormat 'bigwig'  \
	  --binSize 1 --scaleFactor 1.0 --minMappingQuality '1' --samFlagInclude 64 --samFlagExclude 16 --Offset 1

	~/.local/bin/bamCoverage --numberOfProcessors max/2 --bam ${prefix}.proper.q30.bam --outFileName ${prefix}.proper.q30_bottom.bw --outFileFormat 'bigwig'  \
	  --binSize 1 --scaleFactor 1.0 --minMappingQuality '1' --samFlagInclude 80 --Offset 1

done

