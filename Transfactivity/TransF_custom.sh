#!/bin/bash
##Important to keep all of the relative file paths. For some reason Transfactivity absolutely cannot handle spaces in the file paths
##After getting output, the R script annotates it

exp_file=$1
out_dir=${1}_yetfasco_and_custom_Transfactivity

mkdir $out_dir

bin/Transfactivity -expression="$exp_file" -sequence="data/sequence/YeastUpstream.fasta" -psam_list="yetfasco_and_custom.list" -output=$out_dir


Rscript annotate_TF_output.R $out_dir "transfactivity_mvar_coeff.tsv"
Rscript annotate_TF_output.R $out_dir "transfactivity_mvar_pvals.tsv"


