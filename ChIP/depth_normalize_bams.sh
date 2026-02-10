#!/bin/bash

for i in `find . -type f -name "*proper.q30.bam"`
do
	prefix=${i/.bam/}
        echo $prefix
	~/.local/bin/bamCoverage --numberOfProcessors max/2 --bam ${prefix}.bam --outFileName ${prefix}.dnorm.bw --outFileFormat 'bigwig' \
	  --binSize 1 --normalizeUsing CPM -e
done

