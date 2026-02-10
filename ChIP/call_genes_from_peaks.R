


#setwd("~/OneDrive - National Institutes of Health/data/kar4_chip_exo_mr/")


#infile <- "analysis/4_1_Ste12_AF_v_mitotic.txt_win24/4_1_Ste12_AF_v_mitotic.txt_win24_4_1_Ste12_AF_v_mitotic.events.curated.txt"
#infile <- "analysis/1_1_Kar4dep_Ste12_AF.txt_win24/1_1_Kar4dep_Ste12_AF.txt_win24_1_1_Kar4dep_Ste12_AF.events.curated.txt"
infile <- commandArgs(trailingOnly=TRUE)[1]

peaks <- read.delim(file=infile, stringsAsFactors = F)

#optionally filter based on logQ, etc.
# qfilt <- -20
# peaks <- peaks[which(peaks[,grep(pattern = "log2Q", x=colnames(peaks))] <= qfilt),]
#peaks <- peaks[which(is.na(peaks$in_feature)),] #this line would remove any peaks inside of genes/LTRs

#look for instances where a peak should be assigned multiple features, then duplicate it
dcut <- 500 #500 bp cutoff
namevec <- c("genetop", "genebot", "ltrtop", "ltrbot", "in_feature")
namevec2 <- c("genetopdist", "genebotdist", "ltrtopdist", "ltrbotdist", "in_feature_prom_dist")

#could do an apply, but just loop through rows, get nrows now as I will be appending
stoprow <- nrow(peaks)
for (row in 1:stoprow) {
  keep <- which(peaks[row, namevec2] < dcut)
  if (length(keep) > 0) {
    
    #exclude closest
    disc <- which(peaks[row, namevec] == peaks$closest[row])
    keep2 <- setdiff(keep, disc)
    
    #make a new row per remaining, update info
    if (length(keep2) > 0) {
      for (k in keep2) {
        peaks <- rbind(peaks, peaks[row,])
        peaks$closest[nrow(peaks)] <- peaks[row, namevec[k]]
        peaks$closestdist[nrow(peaks)] <- peaks[row, namevec2[k]]
      }
    }
  }
}



#collapse down to just 1 row per feature, keep only a few measures (log2 FC, Q, signal reads, ctrl reads)
#Note, if a feature has multiple peaks associated with it, must decide where to take relevant stats from. Decided to use stats associated with lowest Q-value
keeprows <- c()
likelycats <- c()

#for compatibility with downstream scripts, have to do some column renaming
colnames(peaks)[grep(pattern = "log2Fold", colnames(peaks))] <- "experiment_log2Fold"
colnames(peaks)[grep(pattern = "log2Q", colnames(peaks))] <- "experiment_log2Q"
colnames(peaks)[grep(pattern = "_Sig", colnames(peaks))] <- "experiment_Sig"
colnames(peaks)[grep(pattern = "_Ctrl", colnames(peaks))] <- "experiment_Ctrl"

for (gene in unique(peaks$closest)) {
  generows <- which(peaks$closest==gene)
  keeprow <- generows[which.min(peaks$experiment_log2Q[generows])]
  keeprows <- c(keeprows, keeprow)
  #in keeping same likelycat strategy as before, make LTR if ANY generows have LTR
  likelycats <- c(likelycats, ifelse(any(peaks$likelycat[generows]=="LTR"), "LTR", "gene"))
}

genes <- peaks[keeprows,]

#replace likelycat. comment this out to retain likelycat only for this peak
genes$likelycat <- likelycats 
#genes$likelycat2 <- likelycats #for comparing before/after

#keep same columns as before
genes <- genes[,c("closest", "experiment_log2Fold", "experiment_log2Q", "experiment_Sig", "experiment_Ctrl", 
                    "likelycat", "closestdist")]

genes <- genes[order(genes$experiment_log2Q),]

#could save genes dataframe here, but only has systematic ORF names

#### add systematic ORF name column and "simple" classification for paper ####
#get systematic names for gene list comparison
#load gff annotation
gtf <- read.delim(file="~/OneDrive - National Institutes of Health/data/S288C_reference_genome_R64-4-1_20230830/saccharomyces_cerevisiae_R64-4-1_20230830.gff", 
                  stringsAsFactors = F, skip=21, header = F, quote="")

#has fasta sequence at end marked by ###, don't use
gtf <- gtf[1:(which(gtf[,1] == "###")-1),]

#give column names
colnames(gtf) <- c("chr", "source", "type", "start", "stop", "6", "strand", "introns", "info")

#now create separate columns for various names
gtf$name <- sapply(gtf$info, FUN = function(x){strsplit(x, "Name=")[[1]][2]})
gtf$name <- sapply(gtf$name, FUN = function(x){strsplit(x, ";")[[1]][1]})

gtf$id <- sapply(gtf$info, FUN = function(x){strsplit(x, "ID=")[[1]][2]})
gtf$id <- sapply(gtf$id, FUN = function(x){strsplit(x, ";")[[1]][1]})

gtf$gene <- sapply(gtf$info, FUN = function(x){strsplit(x, "gene=")[[1]][2]})
gtf$gene <- sapply(gtf$gene, FUN = function(x){strsplit(x, ";")[[1]][1]})

#now merge in ORFs for gene names. I built the newgenes annotations so that if there wasn't a gene entry then it uses systematic name already
gtf2 <- gtf[,c("gene", "name")]
gtf2 <- gtf2[which(!(is.na(gtf2$gene))),]

newgenes <- merge(genes, gtf2, by.x="closest", by.y="gene", all.x=T, all.y=F)

#now fill in missing. Could have done this to gtf2 initially
newgenes$name[which(is.na(newgenes$name))] <- newgenes$closest[which(is.na(newgenes$name))]

#now, as my prior analyses only used "simple" genes, add that restriction here
#get orfmap used for earlier analyses
orfmap <- read.delim(file="~/OneDrive - National Institutes of Health/data/yeastract_genemap.txt", stringsAsFactors = F)
#note that when I defined simple genes, I additionally removed mito genes and YIL082W, which is a TY element but was listed as a gene anyway for some reason (maybe the idea was to maintain a single representative TY ORF?)

newgenes$simple <- NA
newgenes$simple[which(newgenes$name %in% orfmap$ORF)] <- 1

newgenes <- newgenes[order(newgenes$simple, newgenes$experiment_log2Q),]

outfile <- paste0(infile, ".genes.txt")
#outfile <- paste0(infile, "_noinfeatures.genes.txt")
write.table(newgenes, file = outfile, sep="\t", row.names = F, quote=F)








