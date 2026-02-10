
#!/bin/bash

chex=~/bin/chexmix.v0.52.public.jar
ginfo=~/builds/chexmix/sacCer3.info
seq=~/builds/bwa/genomes/sacCer3.fa
back=~/builds/chexmix/yeast.back
excludes=~/builds/chexmix/sacCer3_excludes.txt

#bam=36203_Ste12_i5006_MY16679_Ste12-TAP_YPDpH3.5_Alpha-factor_XO.sorted.filtered.proper.q30.bam
#bam=36205_Ste12_i5006_MY16682kar4--kanMX_Ste12-TAP_YPDpH3.5_Alpha-factor_XO_FilteredBAM.proper.q30.bam
#bam2=39270_Ste12_i5006_MY16679_Ste12-TAP_YPDpH3.5_Alpha-factor_XO_FilteredBAM.proper.q30.bam
#bam=39272_Ste12_i5006_MY16682kar4--kanMX_Ste12-TAP_YPDpH3.5_Alpha-factor_XO_FilteredBAM.proper.q30.bam
bam=39270_Ste12_i5006_MY16679_Ste12-TAP_YPDpH3.5_Alpha-factor_XO_FilteredBAM.bam


#for bam in *FilteredBAM.proper.q30.bam
#do
	echo $bam
	out=${bam}_win24

	java -jar $chex --geninfo $ginfo --seq $seq --back $back --expt $bam --format BAM --out $out --threads 8 --verbose --memepath /Users/rogersjv/bin/ --exclude $excludes --mememaxw 24 \
	--round 5 > ${out}.out 2>&1 &
	#--minroc 0.60
	#--betascale 0.05
	#--nomotifs

#done

