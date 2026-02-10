
#!/bin/bash

chex=~/bin/chexmix.v0.52.public.jar
ginfo=~/builds/chexmix/sacCer3.info
seq=~/builds/bwa/genomes/sacCer3.fa
back=~/builds/chexmix/yeast.back
excludes=~/builds/chexmix/sacCer3_excludes.txt

design=1_1_Kar4dep_Ste12_AF.txt

#for design in *txt
#do

	echo $design
	out=${design}_win24

	java -jar $chex --geninfo $ginfo --seq $seq --back $back --design $design --out $out --threads 8 --verbose --memepath /Users/rogersjv/bin/ --exclude $excludes --mememaxw 24 \
	--round 5 > ${out}.out 2>&1 &
#done

