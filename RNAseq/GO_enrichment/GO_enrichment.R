

setwd("~/OneDrive - National Institutes of Health/analysis/kar4_genelist_repeats/")
library(clusterProfiler)
library(org.Sc.sgd.db)


##set up functions

#for retrieving gene names from ensembl
get_gene_names_from_gene_ids <- function(ensemble_ids, annotation_db, 
                                         look_for = 'ENSEMBL', 
                                         fetch = 'GENENAME', return_all = TRUE){
  # Reference organism: Saccharomyces cerevisiae => DATABASE = org.Sc.sgd.db
  symbols <- mapIds(annotation_db, keys = ensemble_ids, column = fetch, 
                    keytype = look_for, multiVals = "first")
  symbols <- symbols[!is.na(symbols)]
  to_name <- ensemble_ids %in% names(symbols)
  ensemble_ids[to_name] <- as.vector(symbols)
  if (return_all){
    return(ensemble_ids)
  }
  else {
    return(ensemble_ids[to_name])
  }
}

# Aux function to convert string of ensembl IDs to string of gene symbol IDs.
get_gene_names <- function(A){
  if(length(A) > 0){
    my_ensembl_ids <- stringr::str_split(string = A, pattern = "/", simplify = T) %>% 
      as.vector(.)
    get_gene_names_from_gene_ids(ensemble_ids = my_ensembl_ids, 
                                 annotation_db = org.Sc.sgd.db, 
                                 look_for = "ORF", fetch = "GENENAME" ) %>%
      stringr::str_flatten(collapse = "/")
  }
  else {
    return("None")
  }
}


## loops over lists of genes, does GO enrichment for all
#load gene lists file
genelists <- read.delim(file="~/OneDrive - National Institutes of Health/paper_figures/alpha_factor_DE/integrated/sumdf_binary_gene_overlap_for_file_S3.txt", stringsAsFactors = F)

#background gene list. Given that all our analyses will use simple + t60-excluded, that will be the universe
uni <- genelists$ORF[which(genelists$simple==1 & genelists$t60_excluded_up==0 & genelists$t60_excluded_down==0)]


#prob here just make a column name vector, then loop over and dump results independently. No real reason to store in a list object.
resnames <- c("wt_up", "wt_down", "kar4_up", "kar4_down", "ste12_up", "ste12_down", "pheromone_induced", 
              "kar4.independent", "kar4.dependent", "kar4D.only")

for (i in 1:length(resnames)) {
  
  #get gene list
  genes <- genelists$ORF[which(genelists[,resnames[i]] == 1)]
  genes <- intersect(genes, uni)
  
  # Overrepresentation analysis for gene list
  go_overrep <- enrichGO(keyType = "ORF",
                            gene = genes,
                            universe      = uni,
                            OrgDb         = org.Sc.sgd.db,
                            ont           = "ALL", #can be "ALL" for all three, or can specify molecular function, cell component, etc.
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = FALSE)
  
  go_overrep.df <- as.data.frame(go_overrep)
  
  #filter out any sets without at least 3 genes
  sizecut <- 3
  go_overrep.df <- go_overrep.df[which(go_overrep.df$Count >= sizecut),]
  
  go_overrep.df$gene_name <- unname(unlist(sapply(go_overrep.df$geneID, function(x){get_gene_names(x)})))
  
  write.table(go_overrep.df, file=paste0("GO_results/go_", resnames[i], ".txt"), row.names = F, quote=F, sep="\t")
  
}




