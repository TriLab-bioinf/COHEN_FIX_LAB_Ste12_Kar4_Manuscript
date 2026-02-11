

library(DESeq2)
options(scipen=999)

setwd("~/OneDrive - National Institutes of Health/analysis/kar4_genelist_repeats/")


#### load data, set up metadata ####

#get gene lists for swapping later
genelist <- read.delim(file="../../data/ORF_gene_map.txt", stringsAsFactors = F)

#get list of simple protein-coding genes for subsetting later. Note there are 7 genes in this list not included in our analysis, I think they were all added to SGD relatively recently
simpgenes <- read.delim(file="../../data/sgd_gene_list.txt", stringsAsFactors = F, header = F)
#append RME2 to simpgenes list specifically for this analysis
simpgenes <- rbind(simpgenes, "RME2")

# Import read counts table into read_counts variable
read_counts.all <- read.table(file = "../../data/raw/read_counts_table_comp8901vs8345.txt",  
                              header = TRUE, # First line contains column headers 
                              row.names = 1, # Column 1 contains row names
                              sep = "\t")   # Column delimiter is a tab (\t)

# Round read counts to the closest integer
read_counts.all <- round(read_counts.all, digits = 0)

# Read metadata.all table
metadata.all <- read.table(file = "../../data/raw/metadata_OC8901vsOC8345.txt", 
                           header = TRUE, 
                           row.names = 1,
                           sep = "\t")

# Sort read_counts.all table so metadata.all and read_counts.all match order (drops a column (length))
if(!unique(colnames(read_counts.all) == rownames(metadata.all))){
  read_counts.all <- read_counts.all[, match(rownames(metadata.all), colnames(read_counts.all))]
}

# Include sample ids in metadata.all as a variable (column)
metadata.all$sample_id <- c("AY31", "AY32" , "AY33" , "AY34" , "AY35" , "AY36" , "AY37" , "AY38" , "AY39" , "AY40" , "AY41" , "AY42", "AY1", "AY2", "AY3", "AY4","AY5", "AY6", "AY7", "AY8", "AY9", "AY10" , "AY11" , "AY12")

# Include total read counts in metadata.all
metadata.all$read_counts <- colSums(read_counts.all, na.rm = TRUE)

# Make group, strain, treatment and genotype columns as factors
metadata.all$treatment <- as.factor(metadata.all$treatment)
metadata.all$strain <- as.factor(metadata.all$strain)
metadata.all$genotype <- as.factor(metadata.all$genotype)
metadata.all$experiment <- as.factor(metadata.all$experiment)


#### could do next sections in one for loop but since there are just a few contrasts doing it manually ####

#store results in this list
res <- vector("list")

#### run DESeq for kar4 t60 vs. wt t60 ####

# Keep only alpha factor 
keep <- metadata.all$treatment=="alpha_factor"
metadata <- metadata.all[keep,]
read_counts <- read_counts.all[,keep]

# drop ste12
keep <- metadata$genotype!="ste12"
metadata <- metadata[keep,]
read_counts <- read_counts[,keep]


##run DESeq
# Add extra metadata column with treatment and genotype data combined
metadata$treat_geno <-  paste(metadata$treatment, metadata$genotype, sep = '_')

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = read_counts,
                              colData = metadata,
                              design = ~ experiment + treat_geno) # Here we assume that all genotypes respond the same to treatment (I think this means to say all genotypes have the same batch effect)
#design = ~ experiment:treat_geno + treat_geno)

# Set reference group for log2FC calculation
#dds$treatment <- relevel(dds$treatment, "none")
dds$genotype <- relevel(dds$genotype, "wt")
dds$treat_geno <- relevel(dds$treat_geno, "alpha_factor_wt") 
dds <- DESeq(dds)
#dds <- DESeq(dds, test="LRT", reduced = ~ experiment)

#### store output
#print(resultsNames(dds))
coefs <- resultsNames(dds)

#coef <- c("treat_geno_alpha_factor_kar4_vs_alpha_factor_wt")
coef <- coefs[length(coefs)]

# Get DESeq2 results. Note that when using coefficient names, you have to convert the coef vector to a list by doing list(c(coef))
DE_results <- results(dds, contrast=list(c(coef)), tidy=T)

# Get shrunk Log2FC
reslfc <- lfcShrink(dds, 
                    contrast = list(c(coef)),
                    type = "ashr")

#can verify p-values are identical
#which(DE_results$padj != reslfc$padj)
#View(cbind(DE_results, reslfc$padj))

#append log2fc to first results
DE_results <- cbind(DE_results, reslfc$log2FoldChange)

#update column names of gene list for later merging. Also merge in gene names here
colnames(DE_results)[ncol(DE_results)] <- "l2fc_shrunk_ashr"
colnames(DE_results)[1] <- "ORF"
DE_results <- merge(genelist, DE_results, by="ORF")
colnames(DE_results) <- paste(colnames(DE_results), coef, sep="_")

res[[1]] <- DE_results
#names(res)[1] <- coef #don't set names as it will append to column names when calling cbind later. Odd... seems like new behvaior

#### run DESeq for wt t60 vs. wt t0 ####

# Keep only WT
keep <- metadata.all$genotype=="wt"
metadata <- metadata.all[keep,]
read_counts <- read_counts.all[,keep]


##run DESeq
# Add extra metadata column with treatment and genotype data combined
metadata$treat_geno <-  paste(metadata$treatment, metadata$genotype, sep = '_')

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = read_counts,
                              colData = metadata,
                              design = ~ experiment + treat_geno) # Here we assume that all genotypes respond the same to treatment (I think this means to say all genotypes have the same batch effect)
#design = ~ experiment:treat_geno + treat_geno)

# Make sure reference groups are "none" and "wt". This means those categories will be the denominator when calculating Log2FC.
#dds$treatment <- relevel(dds$treatment, "none")
dds$genotype <- relevel(dds$genotype, "wt")
dds$treat_geno <- relevel(dds$treat_geno, "none_wt") 
dds <- DESeq(dds)
#dds <- DESeq(dds, test="LRT", reduced = ~ experiment)

#### store output
#print(resultsNames(dds))
coefs <- resultsNames(dds)

#coef <- c("treat_geno_alpha_factor_kar4_vs_alpha_factor_wt")
coef <- coefs[length(coefs)]

# Get DESeq2 results. Note that when using coefficient names, you have to convert the coef vector to a list by doing list(c(coef))
DE_results <- results(dds, contrast=list(c(coef)), tidy=T)

# Get shrunk Log2FC
reslfc <- lfcShrink(dds, 
                    contrast = list(c(coef)),
                    type = "ashr")

#can verify p-values are identical
#which(DE_results$padj != reslfc$padj)
#View(cbind(DE_results, reslfc$padj))

#append log2fc to first results
DE_results <- cbind(DE_results, reslfc$log2FoldChange)

#update column names of gene list for later merging. Also merge in gene names here
colnames(DE_results)[ncol(DE_results)] <- "l2fc_shrunk_ashr"
colnames(DE_results)[1] <- "ORF"
DE_results <- merge(genelist, DE_results, by="ORF")
colnames(DE_results) <- paste(colnames(DE_results), coef, sep="_")

res[[2]] <- DE_results
#names(res)[2] <- coef



#### run DESeq for kar4 t60 vs. kar4 t0 ####

# Keep only kar4
keep <- metadata.all$genotype=="kar4"
metadata <- metadata.all[keep,]
read_counts <- read_counts.all[,keep]


##run DESeq
# Add extra metadata column with treatment and genotype data combined
metadata$treat_geno <-  paste(metadata$treatment, metadata$genotype, sep = '_')

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = read_counts,
                              colData = metadata,
                              design = ~ experiment + treat_geno) # Here we assume that all genotypes respond the same to treatment (I think this means to say all genotypes have the same batch effect)
#design = ~ experiment:treat_geno + treat_geno)

# Make sure reference groups are "none" and "wt". This means those categories will be the denominator when calculating Log2FC.
#dds$treatment <- relevel(dds$treatment, "none")
dds$genotype <- relevel(dds$genotype, "kar4")
dds$treat_geno <- relevel(dds$treat_geno, "none_kar4") 
dds <- DESeq(dds)
#dds <- DESeq(dds, test="LRT", reduced = ~ experiment)

#### store output
#print(resultsNames(dds))
coefs <- resultsNames(dds)

#coef <- c("treat_geno_alpha_factor_kar4_vs_alpha_factor_wt")
coef <- coefs[length(coefs)]

# Get DESeq2 results. Note that when using coefficient names, you have to convert the coef vector to a list by doing list(c(coef))
DE_results <- results(dds, contrast=list(c(coef)), tidy=T)

# Get shrunk Log2FC
reslfc <- lfcShrink(dds, 
                    contrast = list(c(coef)),
                    type = "ashr")

#can verify p-values are identical
#which(DE_results$padj != reslfc$padj)
#View(cbind(DE_results, reslfc$padj))

#append log2fc to first results
DE_results <- cbind(DE_results, reslfc$log2FoldChange)

#update column names of gene list for later merging. Also merge in gene names here
colnames(DE_results)[ncol(DE_results)] <- "l2fc_shrunk_ashr"
colnames(DE_results)[1] <- "ORF"
DE_results <- merge(genelist, DE_results, by="ORF")
colnames(DE_results) <- paste(colnames(DE_results), coef, sep="_")

res[[3]] <- DE_results
#names(res)[3] <- coef

#### run DESeq for ste12 t60 vs. ste12 t0 ####

# Keep only ste12
keep <- metadata.all$genotype=="ste12"
metadata <- metadata.all[keep,]
read_counts <- read_counts.all[,keep]


##run DESeq
# Add extra metadata column with treatment and genotype data combined
metadata$treat_geno <-  paste(metadata$treatment, metadata$genotype, sep = '_')

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = read_counts,
                              colData = metadata,
                              design = ~ experiment + treat_geno) # Here we assume that all genotypes respond the same to treatment (I think this means to say all genotypes have the same batch effect)
#design = ~ experiment:treat_geno + treat_geno)

# Make sure reference groups are "none" and "wt". This means those categories will be the denominator when calculating Log2FC.
#dds$treatment <- relevel(dds$treatment, "none")
dds$genotype <- relevel(dds$genotype, "ste12")
dds$treat_geno <- relevel(dds$treat_geno, "none_ste12") 
dds <- DESeq(dds)
#dds <- DESeq(dds, test="LRT", reduced = ~ experiment)

#### store output
#print(resultsNames(dds))
coefs <- resultsNames(dds)

#coef <- c("treat_geno_alpha_factor_kar4_vs_alpha_factor_wt")
coef <- coefs[length(coefs)]

# Get DESeq2 results. Note that when using coefficient names, you have to convert the coef vector to a list by doing list(c(coef))
DE_results <- results(dds, contrast=list(c(coef)), tidy=T)

# Get shrunk Log2FC
reslfc <- lfcShrink(dds, 
                    contrast = list(c(coef)),
                    type = "ashr")

#can verify p-values are identical
#which(DE_results$padj != reslfc$padj)
#View(cbind(DE_results, reslfc$padj))

#append log2fc to first results
DE_results <- cbind(DE_results, reslfc$log2FoldChange)

#update column names of gene list for later merging. Also merge in gene names here
colnames(DE_results)[ncol(DE_results)] <- "l2fc_shrunk_ashr"
colnames(DE_results)[1] <- "ORF"
DE_results <- merge(genelist, DE_results, by="ORF")
colnames(DE_results) <- paste(colnames(DE_results), coef, sep="_")

res[[4]] <- DE_results
#names(res)[3] <- coef

#### run DESeq for ste12 t0 vs. wt t0 ####

# Keep only no alpha factor 
keep <- metadata.all$treatment=="none"
metadata <- metadata.all[keep,]
read_counts <- read_counts.all[,keep]

# drop kar4
keep <- metadata$genotype!="kar4"
metadata <- metadata[keep,]
read_counts <- read_counts[,keep]


##run DESeq
# Add extra metadata column with treatment and genotype data combined
metadata$treat_geno <-  paste(metadata$treatment, metadata$genotype, sep = '_')

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = read_counts,
                              colData = metadata,
                              design = ~ experiment + treat_geno) # Here we assume that all genotypes respond the same to treatment (I think this means to say all genotypes have the same batch effect)
#design = ~ experiment:treat_geno + treat_geno)

# Make sure reference groups are "none" and "wt". This means those categories will be the denominator when calculating Log2FC.
#dds$treatment <- relevel(dds$treatment, "none")
dds$genotype <- relevel(dds$genotype, "wt")
dds$treat_geno <- relevel(dds$treat_geno, "none_wt") 
dds <- DESeq(dds)
#dds <- DESeq(dds, test="LRT", reduced = ~ experiment)

#### store output
#print(resultsNames(dds))
coefs <- resultsNames(dds)

#coef <- c("treat_geno_alpha_factor_kar4_vs_alpha_factor_wt")
coef <- coefs[length(coefs)]

# Get DESeq2 results. Note that when using coefficient names, you have to convert the coef vector to a list by doing list(c(coef))
DE_results <- results(dds, contrast=list(c(coef)), tidy=T)

# Get shrunk Log2FC
reslfc <- lfcShrink(dds, 
                    contrast = list(c(coef)),
                    type = "ashr")

#can verify p-values are identical
#which(DE_results$padj != reslfc$padj)
#View(cbind(DE_results, reslfc$padj))

#append log2fc to first results
DE_results <- cbind(DE_results, reslfc$log2FoldChange)

#update column names of gene list for later merging. Also merge in gene names here
colnames(DE_results)[ncol(DE_results)] <- "l2fc_shrunk_ashr"
colnames(DE_results)[1] <- "ORF"
DE_results <- merge(genelist, DE_results, by="ORF")
colnames(DE_results) <- paste(colnames(DE_results), coef, sep="_")

res[[5]] <- DE_results
#names(res)[1] <- coef #don't set names as it will append to column names when calling cbind later. Odd... seems like new behvaior

#### run DESeq for kar4 t0 vs. wt t0 ####

# Keep only no alpha factor 
keep <- metadata.all$treatment=="none"
metadata <- metadata.all[keep,]
read_counts <- read_counts.all[,keep]

# drop ste12
keep <- metadata$genotype!="ste12"
metadata <- metadata[keep,]
read_counts <- read_counts[,keep]


##run DESeq
# Add extra metadata column with treatment and genotype data combined
metadata$treat_geno <-  paste(metadata$treatment, metadata$genotype, sep = '_')

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = read_counts,
                              colData = metadata,
                              design = ~ experiment + treat_geno) # Here we assume that all genotypes respond the same to treatment (I think this means to say all genotypes have the same batch effect)
#design = ~ experiment:treat_geno + treat_geno)

# Make sure reference groups are "none" and "wt". This means those categories will be the denominator when calculating Log2FC.
#dds$treatment <- relevel(dds$treatment, "none")
dds$genotype <- relevel(dds$genotype, "wt")
dds$treat_geno <- relevel(dds$treat_geno, "none_wt") 
dds <- DESeq(dds)
#dds <- DESeq(dds, test="LRT", reduced = ~ experiment)

#### store output
#print(resultsNames(dds))
coefs <- resultsNames(dds)

#coef <- c("treat_geno_alpha_factor_kar4_vs_alpha_factor_wt")
coef <- coefs[length(coefs)]

# Get DESeq2 results. Note that when using coefficient names, you have to convert the coef vector to a list by doing list(c(coef))
DE_results <- results(dds, contrast=list(c(coef)), tidy=T)

# Get shrunk Log2FC
reslfc <- lfcShrink(dds, 
                    contrast = list(c(coef)),
                    type = "ashr")

#can verify p-values are identical
#which(DE_results$padj != reslfc$padj)
#View(cbind(DE_results, reslfc$padj))

#append log2fc to first results
DE_results <- cbind(DE_results, reslfc$log2FoldChange)

#update column names of gene list for later merging. Also merge in gene names here
colnames(DE_results)[ncol(DE_results)] <- "l2fc_shrunk_ashr"
colnames(DE_results)[1] <- "ORF"
DE_results <- merge(genelist, DE_results, by="ORF")
colnames(DE_results) <- paste(colnames(DE_results), coef, sep="_")

res[[6]] <- DE_results
#names(res)[1] <- coef #don't set names as it will append to column names when calling cbind later. Odd... seems like new behvaior

#### run DESeq for ste12 t60 vs. wt t60 ####

# Keep only alpha factor 
keep <- metadata.all$treatment=="alpha_factor"
metadata <- metadata.all[keep,]
read_counts <- read_counts.all[,keep]

# drop kar4
keep <- metadata$genotype!="kar4"
metadata <- metadata[keep,]
read_counts <- read_counts[,keep]


##run DESeq
# Add extra metadata column with treatment and genotype data combined
metadata$treat_geno <-  paste(metadata$treatment, metadata$genotype, sep = '_')

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = read_counts,
                              colData = metadata,
                              design = ~ experiment + treat_geno) # Here we assume that all genotypes respond the same to treatment (I think this means to say all genotypes have the same batch effect)
#design = ~ experiment:treat_geno + treat_geno)

# Make sure reference groups are "none" and "wt". This means those categories will be the denominator when calculating Log2FC.
#dds$treatment <- relevel(dds$treatment, "none")
dds$genotype <- relevel(dds$genotype, "wt")
dds$treat_geno <- relevel(dds$treat_geno, "alpha_factor_wt") 
dds <- DESeq(dds)
#dds <- DESeq(dds, test="LRT", reduced = ~ experiment)

#### store output
#print(resultsNames(dds))
coefs <- resultsNames(dds)

#coef <- c("treat_geno_alpha_factor_kar4_vs_alpha_factor_wt")
coef <- coefs[length(coefs)]

# Get DESeq2 results. Note that when using coefficient names, you have to convert the coef vector to a list by doing list(c(coef))
DE_results <- results(dds, contrast=list(c(coef)), tidy=T)

# Get shrunk Log2FC
reslfc <- lfcShrink(dds, 
                    contrast = list(c(coef)),
                    type = "ashr")

#can verify p-values are identical
#which(DE_results$padj != reslfc$padj)
#View(cbind(DE_results, reslfc$padj))

#append log2fc to first results
DE_results <- cbind(DE_results, reslfc$log2FoldChange)

#update column names of gene list for later merging. Also merge in gene names here
colnames(DE_results)[ncol(DE_results)] <- "l2fc_shrunk_ashr"
colnames(DE_results)[1] <- "ORF"
DE_results <- merge(genelist, DE_results, by="ORF")
colnames(DE_results) <- paste(colnames(DE_results), coef, sep="_")

res[[7]] <- DE_results
#names(res)[1] <- coef #don't set names as it will append to column names when calling cbind later. Odd... seems like new behvaior


#### aggregate data and compare to prior gene lists, save for motif analyses and heatmaps ####

#aggregate
resdf <- do.call(cbind, res)

#save results
write.table(resdf, file="all_contrasts.txt", row.names = F, sep="\t", quote=F)

#subset on simple (protein coding gene from SGD) gene list
resdf <- resdf[which(resdf$ORF_treat_geno_alpha_factor_kar4_vs_alpha_factor_wt %in% simpgenes[,1]),]

#save simple results
write.table(resdf, file="all_contrasts_simple.txt", row.names = F, sep="\t", quote=F)



