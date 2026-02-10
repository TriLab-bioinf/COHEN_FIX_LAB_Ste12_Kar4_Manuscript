**File list and explanation:**

**TransF_custom.sh:** Provides code to run Transfactivity and annotate output. Takes as input an expression file (log2-FC matrix, one header row and first column is ORF gene names (must be systematic ORF names)), a FASTA file containing gene promoter sequences (provided with REDUCE suite), and a file listing all the PSAMs (motifs) to use (provided here as yetfasco_and_custom.list; note individual motifs can be commented out using #). Assumes the script is being run from within the REDUCE suite folder (http://reducesuite.bussemakerlab.org/).

**annotate_TF_output.R:** Provides code to annotate the Transfactivity output (primarily converts PSAM IDs to common gene names for human-readability, and converts p-values to signed -log10 p-values). Should not need to be modified, gets called by the TransF_custom.sh script.

**convert2psam.sh:** Provides example of how to run Convert2PSAM for converting .pfm files (for example, those used by YetFaSCo (http://yetfasco.ccbr.utoronto.ca/1.02/Downloads/Expert_PFMs.zip) to PSAM (.xml) format. Examples of custom .pfm and converted files are also given (H-T 4 and T-T 3 motifs derived from MEME). Note that as Transfactivity was run in multivariate mode, results will not be the same unless all of the yetfasco motifs are also included in the analysis. 


**NOTE:** These scripts provide the explicit commands used to perform the analyses described in the paper. However, they are not intended to be functional for all users without modification. Users are expected to install the required programs on their own machine (or computing cluster), then update the files to point at their own filepaths.
