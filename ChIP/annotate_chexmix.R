


#setwd("~/OneDrive - National Institutes of Health/data/kar4_chip_exo_mr/")

ncrna <- 0 #if set to 1, will include hernan's ncRNA transcripts as possible genes


#### load data ####

#load binding events list
#infile <- "36203_Ste12_i5006_MY16679_Ste12-TAP_YPDpH3.5_Alpha-factor_XO_FilteredBAM.proper.q30.bam_win24/36203_Ste12_i5006_MY16679_Ste12-TAP_YPDpH3.5_Alpha-factor_XO_FilteredBAM.proper.q30.bam_win24_experiment.events"
#infile <- "analysis/4_1_Ste12_AF_v_mitotic.txt_win24/4_1_Ste12_AF_v_mitotic.txt_win24_4_1_Ste12_AF_v_mitotic.events"
#infile <- "analysis/1_1_Kar4dep_Ste12_AF.txt_win24/1_1_Kar4dep_Ste12_AF.txt_win24.all.events.table"

#use this to call from command line, assumes just passing a single filename
infile <- commandArgs(trailingOnly=TRUE)[1]

#peaks <- read.delim(file=infile, skip=6, header = T) #6 if no reps specified file used
# The # lines to skip varies, so search for it dynamically rather than specify
#It uses max width of first 5 rows to determine column number, I think
#Skip 6 rows as there will always be at least that many rows to skip, but sometimes more
peaks <- read.delim(file=infile, header = F, skip=6)
skipind <- 5 + which(peaks[,1]=="#Point")

#now reload with headers
peaks <- read.delim(file=infile, skip=skipind, header = T)

if (nrow(peaks) > 0 & colnames(peaks)[1]=="X.Point") {
  #separate chr and site info
  peaks$chr <- sapply(peaks$X.Point, FUN = function(x){strsplit(x, ":")[[1]][1]})
  peaks$site <- as.numeric(sapply(peaks$X.Point, FUN = function(x){strsplit(x, ":")[[1]][2]}))
  
  
  #add columns to be filled in later
  peaks$genetop <- NA #closest ORF with correct promoter orientation for top strand
  peaks$genetopdist <- NA
  peaks$genebot <- NA #closest ORF with correct promoter orientation for bottom strand
  peaks$genebotdist <- NA
  peaks$ltrtop <- NA #closest LTR. True TY elements already captured in gene category
  peaks$ltrtopdist <- NA
  peaks$ltrbot <- NA #closest LTR. True TY elements already captured in gene category
  peaks$ltrbotdist <- NA
  peaks$closest <- NA #aggregates closest of the 4 categories for ease of viewing
  peaks$closestdist <- NA
  peaks$in_feature <- NA #for when peak is in a gene/feature
  peaks$in_feature_prom_dist <- NA
  peaks$likelycat <- NA #for marking when LTR is closest or in_feature. Not perfect, but want to throw these away for fasta generation
  
  
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
  
  
  
  #optionally merge in Hernan's ncRNA transcript list here
  if (ncrna == 1) {
    ncdf <- read.delim(file="~/OneDrive - National Institutes of Health/data/DE_af_vs_none_ALL.txt", stringsAsFactors = F)
    #discard bottom 2 rows
    ncdf <- ncdf[1:684,]
    
    #for some reason end3 is a string, not numeric. fix
    ncdf$end3 <- as.numeric(ncdf$end3)
    
    #make it look like gtf for merging
    ncgtf <- data.frame(chr=paste0("chr", ncdf$chr),
                        source="StringTie",
                        type="gene",
                        start=ncdf$end5, 
                        stop=ncdf$end3,
                        '6'=".",
                        strand=ncdf$strand,
                        introns=".",
                        info="",
                        name=ncdf$ID,
                        id=ncdf$ID,
                        gene=ncdf$ID,
                        check.names = FALSE)
    
    #rename chrMito to chrM
    ncgtf$chr[which(ncgtf$chr=="chrMito")] <- "chrM"
    
    gtf <- rbind(gtf, ncgtf)
  }
  
  
  #### curate events table ####
  
  for (row in 1:nrow(peaks)) {
    #first get top/bottom candidate gene promoters
    #get all genes on top strand of relevant chromosome and in front of peak site
    topgenes <- which(gtf$type=="gene" & gtf$chr==peaks$chr[row] & gtf$strand=="+" & gtf$start >= peaks$site[row])
    if (length(topgenes) > 0) {
      dists <- gtf$start[topgenes] - peaks$site[row]
      keepgene <- which.min(dists) #they should already be ordered by starting position, but I'm writing the code without making any assumptions. They could be randomized even.
      peaks$genetopdist[row] <- dists[keepgene]
      peaks$genetop[row] <- gtf$gene[topgenes[keepgene]]
      if (is.na(peaks$genetop[row])) {
        peaks$genetop[row] <- gtf$name[topgenes[keepgene]]
      }
    }
    
    #get all genes on bottom strand of relevant chromosome and in front of peak site
    botgenes <- which(gtf$type=="gene" & gtf$chr==peaks$chr[row] & gtf$strand=="-" & gtf$stop <= peaks$site[row])
    if (length(botgenes) > 0) {
      dists <- peaks$site[row] - gtf$stop[botgenes]
      keepgene <- which.min(dists) #they should already be ordered by starting position, but I'm writing the code without making any assumptions. They could be randomized even.
      peaks$genebotdist[row] <- dists[keepgene]
      peaks$genebot[row] <- gtf$gene[botgenes[keepgene]]
      if (is.na(peaks$genebot[row])) {
        peaks$genebot[row] <- gtf$name[botgenes[keepgene]]
      }
    }
    
    
    #next get top/bottom candidate LTR promoters
    #get all ltrs on top strand of relevant chromosome and in front of peak site
    topltrs <- which(gtf$type=="long_terminal_repeat" & gtf$chr==peaks$chr[row] & gtf$strand=="+" & gtf$start >= peaks$site[row])
    if (length(topltrs) > 0) {
      dists <- gtf$start[topltrs] - peaks$site[row]
      keepltr <- which.min(dists) #they should already be ordered by starting position, but I'm writing the code without making any assumptions. They could be randomized even.
      peaks$ltrtopdist[row] <- dists[keepltr]
      peaks$ltrtop[row] <- gtf$name[topltrs[keepltr]]
    }
    
    #get all ltrs on bottom strand of relevant chromosome and in front of peak site
    botltrs <- which(gtf$type=="long_terminal_repeat" & gtf$chr==peaks$chr[row] & gtf$strand=="-" & gtf$stop <= peaks$site[row])
    if (length(botltrs) > 0) {
      dists <- peaks$site[row] - gtf$stop[botltrs]
      keepltr <- which.min(dists) #they should already be ordered by starting position, but I'm writing the code without making any assumptions. They could be randomized even.
      peaks$ltrbotdist[row] <- dists[keepltr]
      peaks$ltrbot[row] <- gtf$name[botltrs[keepltr]]
    }
    
    #now get minimum of all 4 categories.
    #use a column names vector to avoid using column indices
    namevec <- c("genetop", "genebot", "ltrtop", "ltrbot")
    namevec2 <- c("genetopdist", "genebotdist", "ltrtopdist", "ltrbotdist")
    
    keep <- which.min(peaks[row, namevec2])
    
    peaks$closest[row] <- peaks[row, namevec[keep]]
    peaks$closestdist[row] <- peaks[row, namevec2[keep]]
    if (keep %in% c(3,4)) {
      peaks$likelycat[row] <- "LTR"
    } else {
      peaks$likelycat[row] <- "gene"
    }
    
    
    #now test if it's inside a feature
    ingenes <- which((gtf$type=="gene" | gtf$type=="long_terminal_repeat") & gtf$chr==peaks$chr[row] & gtf$start <= peaks$site[row] & gtf$stop >= peaks$site[row])
    if (length(ingenes) > 0) {

      #theoretically could be in multiple features. Use the largest feature in this case
      sizes <- gtf$stop[ingenes] - gtf$start[ingenes]
      keepgene <- which.max(sizes)
      peaks$in_feature[row] <- gtf$gene[ingenes[keepgene]]
      if (is.na(peaks$in_feature[row])) {
        peaks$in_feature[row] <- gtf$name[ingenes[keepgene]]
      }
      
      #get distance from feature start, useful in gene calling script
      if (gtf$strand[ingenes[keepgene]] == "+") {
        peaks$in_feature_prom_dist[row] <- peaks$site[row] - gtf$start[ingenes[keepgene]]
      } else if (gtf$strand[ingenes[keepgene]] == "-") {
        peaks$in_feature_prom_dist[row] <- gtf$stop[ingenes[keepgene]] - peaks$site[row]
      }
      
      if (gtf$type[ingenes[keepgene]] == "long_terminal_repeat") {
        peaks$likelycat[row] <- "LTR"
      }
    }
  }

  if (ncrna==1) {
    outfile <- paste0(infile, "_ncrna.curated.txt")
  } else {
    outfile <- paste0(infile, ".curated.txt")
  }
  
  
  
  write.table(peaks, file = outfile, sep="\t", row.names = F, quote=F)
}

