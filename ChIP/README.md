
**File list and explanation:**

**00_alignreads.sh:** Provides code for initial alignment and filtering of reads, takes as input raw fastq.gz read files and outputs a BAM.

**01_clean_count_bam.sh:** Provides code for using samtools to further clean BAM files and keep only reads with proper pair flag and MQ 30. Also contains optional code for creating tag count bigwigs for investigating tag count distributions. 

**depth_normalize_bams.sh:** Provides code for creating depth-normalized coverage files (bigwig format) for visually comparing relative occupancy in IGV. 

**run_multibigwigsummary.sh:** Provides code for counting tags (from ChIP-exo) over specified regions (from a BED file). Uses the tag count bigwigs from 01_clean_count_bam.sh as input. Note that the raw output file provides counts/length (I.e., the average tag count over a region). To get absolute tag counts, multiply the output by each region length in the BED file.

**02_chexmix.sh:** Provides code to run chexmix on individual bam files (no background control). All files (genomic background nucleotide frequencies, for example) can be found on the ChExMix GitHub, except the region exclusion file which is provided here (in misc_files).

**chexmix_fulldesign.sh:** Provides code to run chexmix using a design file (pools replicates and specifies background controls to perform differential contrasts against). An example design file is provided at design_files/1_1_Kar4dep_Ste12_AF.txt. Other design files used in this study are also included.

**annotate_chexmix.R:** Curates chexmix output to associate peaks with likely gene promoters. Takes as input a .events file from chexmix and outputs a curated .events file. Has the option to include intergenic transcripts in gene searching, use the file provided at misc_files/DE_af_vs_none_ALL.txt.

**call_genes_from_peaks.R:** Converts curated peaks files to a gene-centric list (one row per gene, the same peak can be associated with multiple genes, and if multiple peaks are associated with one gene, it gets reduced down to statistics associated with the strongest peak). Uses a file containing a list of “simple” genes for classification which can be found in misc_files/yeastract_genemap.txt. Note this can also easily be replaced by using the “simple” gene list found in Supplementary File S3.

**convert_peaks_to_fastas.R:** Takes a curated peaks file as input and outputs a DNA-sequence (.fasta format) corresponding to the peak region, which was used as input to the MEME and FIMO online tools. Also outputs a modified curated peaks file for merging with FIMO analyses. Note that for FIMO we also provide the motifs file used (misc_files/ste12_dimotifs_full.txt).

**compare_peaks_to_motifs.R:** For curating FIMO output (# of mismatches, etc.). Takes as input a curated merged peaks file (output from convert_peaks_to_fastas.R), as well as FIMO outputs (examples of these three input files are given in the misc_files folder (files starting with 1_1_Kar4dep…). Note you will have to unzip the fimo.tsv file first). Outputs a file which is identical to the curated merged peaks input file, but with added info related to the best di-motif for that peak region from FIMO.


**NOTE:** These scripts provide the explicit commands used to perform the analyses described in the paper. However, they are not intended to be functional for all users without modification. Users are expected to install the required programs on their own machine (or computing cluster), then update the files to point at their own filepaths.



