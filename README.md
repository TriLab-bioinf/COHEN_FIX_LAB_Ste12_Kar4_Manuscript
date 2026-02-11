# Kar4 acts as a Ste12 regulator in Saccharomyces cerevisiae, promoting Ste12 binding to a specific DNA motif genome-wide

## Jason V. Rogers<sup>1*</sup>, Amanda Yeo<sup>1*</sup>, Val Meleshkevich<sup>1</sup>, Hernan Lorenzi<sup>1**</sup>, Mark D. Rose<sup>2**</sup>, Orna Cohen-Fix<sup>1**</sup>

_1: The Laboratory of Biochemistry and Genetics, National Institute of Diabetes and Digestive and Kidney Disease, National Institutes of Health, Bethesda MD, 20892, USA_

_2: Department of Biology, Georgetown University, Washington DC, 20057 USA_

*Co-first author, **Co-Corresponding authors


[<img width="191" height="20" alt="image" src="https://github.com/user-attachments/assets/d87e4837-da00-42c7-96c8-f77f830b22a8" />](https://doi.org/10.5281/zenodo.17703432)

Folder [ChIP/](ChIP) contains scripts and helper files for ChIP-exo processing and analysis, including read alignment and preprocessing (e.g. `00_alignreads.sh`, `01_clean_count_bam.sh`), downstream analysis with ChExMix (`02_chexmix.sh`), peak annotation and motif comparison scripts.

Folder [RNAseq/](RNAseq) contains pipelines and analysis for RNA-seq data, including read processing (`READ_processing/`), exploratory analyses usinf heatmaps and PCA (`Exploratory_analysis/`), differential expression analysis (`Differential_expression/` with `DESeq.R`), and GO enrichment utilities (`GO_enrichment/`) together with the necessary data and metadata.

Folder [Transfactivity/](Transfactivity) contains tools and data for transcription factor activity analysis, including scripts to annotate TF outputs and convert motif formats (`annotate_TF_output.R`, `convert2psam.sh`, `TransF_custom.sh`), PSAM libraries in `data/PSAMs/`, and helper files such as `SGD_features.tab` and `yetfasco_and_custom.list`.

Folder [Intergenic_transcripts/STRINGTIE_intergenic/](Intergenic_transcripts/STRINGTIE_intergenic) contains pipeline to predict candidate novel intergenic transcripts in _Saccharomyces cerevisiae_ from RNAseq data.

Folder [Intergenic_transcripts/Differencial_expression_intergenic/](Intergenic_transcripts/Differencial_expression_intergenic) contains R markdown script `intergenic_DE_analysis.Rmd` and accessory files for differential expression analysis of predicted intergenic transcripts.

**NOTE:** These scripts provide the explicit commands used to perform the analyses described in the paper (i.e., to convert the raw data to the processed data found in the Supplementary Files). However, they are not intended to be functional for all users without modification. Users are expected to install the required programs on their own machine (or computing cluster), then update the files to point at their own filepaths. Raw RNA-seq read files and the gene-count matrix can be found at GEO accession ID GSE306669. Raw ChIP-exo read files and final depth-normalized read alignments for use with IGV can be found at GEO accession ID GSE306667.
