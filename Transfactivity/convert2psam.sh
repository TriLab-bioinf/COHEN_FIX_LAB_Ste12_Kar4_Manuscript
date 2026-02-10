#!/bin/sh

#for file in `ls data/yetfasco/ALIGNED_ENOLOGO_FORMAT_PFMS/`
for file in `ls data/custom/`
do
	#bin/Convert2PSAM -source=pf -inp=data/yetfasco/ALIGNED_ENOLOGO_FORMAT_PFMS/${file} -psam=data/PSAMs/yetfasco/${file}.xml
	bin/Convert2PSAM -source=pf -inp=data/custom/${file} -psam=data/PSAMs/custom/${file}.xml
done
