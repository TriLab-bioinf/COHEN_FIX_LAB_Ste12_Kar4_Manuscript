
#convert chexmix output peaks to fasta files for use with MEME etc.


#setwd("~/OneDrive - National Institutes of Health/data/kar4_chip_exo_mr/")

#biostrings functions: https://www.bioconductor.org/packages/2.13/bioc/vignettes/Biostrings/inst/doc/BiostringsQuickOverview.pdf
suppressPackageStartupMessages(library(Biostrings))


#### load and prep data ####
#load genome file
#genome <- readDNAStringSet("~/OneDrive - National Institutes of Health/data/S288C_reference_genome_R64-4-1_20230830/S288C_reference_sequence_R64-4-1_20230830.fsa")
genome <- readDNAStringSet(("~/igv/genomes/seq/sacCer3.fa"))


#infile <- "analysis/1_1_Kar4dep_Ste12_AF.txt_win24/1_1_Kar4dep_Ste12_AF.txt_win24_1_1_Kar4dep_Ste12_AF.events.curated.txt"
#use this to call from command line, assumes just passing a single filename
infile <- commandArgs(trailingOnly=TRUE)[1]

#load curated peaks. Script is hard-designed to work with these files now
peaks <- read.delim(file=infile, stringsAsFactors = F)

#keep only non-LTR sequences
peaks <- peaks[which(peaks$likelycat!="LTR"),]


#Examples of optional filtering steps
# #special filtering for all.events file
# peaks <- peaks[which(peaks$X1_1_Kar4dep_Ste12_AF_log2Fold <= 1 & 
#                        peaks$X1_1_Kar4dep_Ste12_AF_log2Fold >= -1 & 
#                        peaks$X1_1_Kar4dep_Ste12_AF_Sig >= 100 &
#                        peaks$X1_1_Kar4dep_Ste12_AF_log2Q >= -5),]

#optionally increase Q or log2FC cutoff
#peaks <- peaks[which(peaks$X1_1_Kar4dep_Ste12_AF_Sig >= 50),]
#peaks <- peaks[which(peaks$X1_1_Kar4dep_Ste12_AF_log2Fold >= 3.15),]
#peaks <- peaks[which(peaks$X1_1_Kar4dep_Ste12_AF_log2Q < -35),]

#optionally take top n peaks. They usually come pre-sorted by peak height.
#peaks <- peaks[1:200,]


#### extract promoter sequences, save  ####

#set how many bp to use (including peak will be bpup + bpdown + 1)
bpup <- as.numeric(commandArgs(trailingOnly=TRUE)[2]) #pass size value from command line (after filename)
#bpup <- 60 #alternatively, comment above line and set the value here
bpdown <- bpup

if (nrow(peaks) > 0) {
  #sort for peak merging later
  peaks <- peaks[order(peaks$chr, peaks$site),]
  
  ##to handle peaks that are too close to each other, I first calculate and save the start/stop sites to extract for fasta file, then merge those that are overlapping, ending up with a new dataframe for fasta generation
  peaks$start <- NA
  peaks$stop <- NA
  peaks$newID <- NA
  peaks$newname <- NA
  
  for (row in 1:nrow(peaks)) {
    start <- peaks$site[row]-bpdown
    stop <- peaks$site[row]+bpup
    
    if (start <= 0) {
      start <- 1
    }
    
    genlen <- length(genome[[peaks$chr[row]]])
    if (stop > genlen) {
      stop <- genlen
    }
    
    #store assignments
    peaks$start[row] <- start
    peaks$stop[row] <- stop
    if (row == 1) {
      peaks$newID[row] <- 1
    } else {
      peaks$newID[row] <- 1 + peaks$newID[row-1]
    }
    peaks$newname[row] <- peaks$X.Point[row]
    
    
    #now check for close peaks to merge. To disable merging just comment out this chunk
    #for now, just merges if they actually overlap. However, might be more sensible to merge if they overlap by a certain amount, like 7 bp or whatever expected motif size is
      #in looking at 4_1 file example, only once is overlap less than 7 bp
    if (row > 1) {
      if ((start < peaks$stop[row-1]) & peaks$chr[row]==peaks$chr[row-1]) {
        #print(row)
        #print(peaks$stop[row-1] - start) #print gap size
        peaks$newID[row] <- peaks$newID[row-1]
        peaks$newname[row] <- paste(peaks$newname[row-1], peaks$newname[row], sep="_")
      }
    }
  }

  #aggregate necessary info. This merges overlapping peaks
  peaksb <- aggregate(start ~ newID, peaks, min)
  peaksc <- aggregate(stop ~ newID, peaks, max)
  peaksd <- aggregate(cbind(chr, newname) ~ newID, peaks, FUN = function(x){x[which.max(nchar(x))]}) #just carries over chr indication
  peaks2 <- merge(peaksb, peaksc, by="newID")
  peaks2 <- merge(peaks2, peaksd, by="newID")
  
  #Carry over some extra info from peaks to peaks2 and save for a downstream script for curating FIMO results
   
  #want to save log2Fold and Q info. Generalized this for any contrast by searching column names
  cname1 <- colnames(peaks)[grep(pattern = "_log2Fold", x=colnames(peaks))[1]]
  peaks2$tmpcolname <- NA #not sure how to assign column name as variable here, so just rename later
  colnames(peaks2)[ncol(peaks2)] <- cname1
  
  cname2 <- colnames(peaks)[grep(pattern = "_log2Q", x=colnames(peaks))[1]]
  peaks2$tmpcolname2 <- NA #not sure how to assign column name as variable here, so just rename later. don't technically need to give this column a new name, could just keep re-using tmpcolname
  colnames(peaks2)[ncol(peaks2)] <- cname2
  
  cname3 <- colnames(peaks)[grep(pattern = "_Sig", x=colnames(peaks))[1]]
  peaks2$tmpcolname3 <- NA #not sure how to assign column name as variable here, so just rename later
  colnames(peaks2)[ncol(peaks2)] <- cname3
  
  cname4 <- colnames(peaks)[grep(pattern = "_Ctrl", x=colnames(peaks))[1]]
  peaks2$tmpcolname4 <- NA #not sure how to assign column name as variable here, so just rename later
  colnames(peaks2)[ncol(peaks2)] <- cname4
  
  for (row in 1:nrow(peaks2)) {
    userows <- which(peaks$newID == peaks2$newID[row])
    userow <- userows[which.min(peaks[userows,cname2])] #gets stats from lowest log2Q peak when merged
    peaks2[row, cname1] <- peaks[userow, cname1]
    peaks2[row, cname2] <- peaks[userow, cname2]
    peaks2[row, cname3] <- peaks[userow, cname3]
    peaks2[row, cname4] <- peaks[userow, cname4]
    peaks2[row, "closest"] <- peaks[userow, "closest"]
    peaks2[row, "closestdist"] <- peaks[userow, "closestdist"]
    peaks2[row, "in_feature"] <- peaks[userow, "in_feature"]
    peaks2[row, "in_feature_prom_dist"] <- peaks[userow, "in_feature_prom_dist"]
  }
  
  #save peaks2
  bpstring <- as.character(bpup+bpdown+1)
  write.table(peaks2, file = paste0(infile, "_dupsmerged_", bpstring,  "bp.txt"), sep="\t", quote=F, row.names = F)
  
  ##now extract sequences for fasta file
  peakseqs <- DNAStringSet()

  #could build a chr index lookup table here, but since I used same genome file for chexmix, they should be identical names
  #since I'm indexing by name, doesn't matter what order they're in, but do note they are in alphabetical not numerical order in genome file
  for (row in 1:nrow(peaks2)) {
    peakseqs[[row]] <- genome[[peaks2$chr[row]]][peaks2$start[row]:peaks2$stop[row]]
  }

  
  #give names, prob just peak names
  names(peakseqs) <- peaks2$newname
  
  #turns out MEME (at least online tool) rejects name strings > 50 characters
  #truncate here
  longnames <- which(nchar(names(peakseqs)) > 50)
  if (length(longnames) > 0) {
    for (name in longnames) {
      names(peakseqs)[name] <- substr(names(peakseqs)[name], 1, 50)
    }
  }
  
  
  #save fasta
  #bpstring <- as.character(bpup+bpdown+1)
  writeXStringSet(peakseqs, filepath = paste0(infile, "_dupsmerged_", bpstring,  "bp.fasta"), format="fasta")
  
  #NOTE, if sending to FIMO there is a bug where it truncates sequence names containing a ":". Thus, rename here
  names(peakseqs) <- gsub(pattern = ":", replacement = "_", x=names(peakseqs))
  writeXStringSet(peakseqs, filepath = paste0(infile, "_dupsmerged_", bpstring,  "bp_forFIMO.fasta"), format="fasta")
}














