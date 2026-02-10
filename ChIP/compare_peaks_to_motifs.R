

setwd("~/OneDrive - National Institutes of Health/data/kar4_chip_exo_mr/")
library(Biostrings)
#library(ggplot2)
#theme_set(theme_bw())


#load a curated peaks file after merging from the convert_peaks_to_fasta.R script
peaks2 <- read.delim(file="analysis/1_1_Kar4dep_Ste12_AF.txt_win24/1_1_Kar4dep_Ste12_AF.txt_win24_1_1_Kar4dep_Ste12_AF.events.curated.txt_dupsmerged_241bp.txt", stringsAsFactors = F)


#### look at fimo output

#start with the "best" file 
fimo <- read.delim(file="analysis/1_1_Kar4dep_Ste12_AF.txt_win24/1_1_Kar4dep_Ste12_AF.txt_win24_1_1_Kar4dep_Ste12_AF.events.curated.txt_dupsmerged_241bp_nc.fasta_fimo_best_site.narrowPeak", stringsAsFactors = F, header = F)

#add annotations
colnames(fimo) <- c("sequence_name", "start", "stop", "motif_id", "noidea", "strand", "score", "p.value", "q.value", "noidea.1")
fimo$orientation <- sapply(fimo$motif_id, FUN = function(x) {strsplit(x, "_")[[1]][1]})
fimo$m_length <- as.numeric(sapply(fimo$motif_id, FUN = function(x) {strsplit(x, "_")[[1]][2]}))

#optionally exclude a motif type here
#fimo <- fimo[which(fimo$motif_id != "H-T_4_nt"),]
#fimo <- fimo[which(fimo$motif_id == "H-T_4_nt"),]
#fimo <- fimo[which(fimo$motif_id == "T-T_3_nt"),]

#keep only best per peak
fimo$best <- 0
seqs <- unique(fimo$sequence_name)
for (seq in seqs) {
  rows <- which(fimo$sequence_name == seq)
  #should already be ordered, so should always end up being the first row, but just do it without making this assumption anyway
  fimo$best[rows[which.min(fimo$q.value[rows])]] <- 1
}

fimo2 <- fimo[which(fimo$best==1),]

# #make a histogram of results
# ggplot(fimo2, aes(x=m_length, fill=orientation)) +
#   geom_histogram(binwidth=1) +
#   facet_wrap(~ orientation)

#View(ftable(fimo2$motif_id))


#calculate # mismatches for each motif
#have to retrieve full fimo output to get exact matched sequence, then subset
motifs <- read.delim(file="analysis/1_1_Kar4dep_Ste12_AF.txt_win24/1_1_Kar4dep_Ste12_AF.txt_win24_1_1_Kar4dep_Ste12_AF.events.curated.txt_dupsmerged_241bp_nc.fasta_fimo.tsv", stringsAsFactors = F)


fimo2$fullID <- paste(fimo2$sequence_name, fimo2$motif_id, fimo2$start, fimo2$stop, fimo2$strand, sep="_")
motifs$fullID <- paste(motifs$sequence_name, motifs$motif_id, motifs$start, (motifs$stop+1), motifs$strand, sep="_") #note that the stop sites are defined differently between files... so I have to add 1

fimo2 <- merge(fimo2, motifs[,c("fullID", "matched_sequence")], by="fullID")

#extract motifs found
fimo2$motif1 <- substr(fimo2$matched_sequence, 1, 7)
fimo2$motif2 <- substr(fimo2$matched_sequence, 8+fimo2$m_length, 14+fimo2$m_length)

#RC motifs for T-T and H-H matches
for (row in 1:nrow(fimo2)) {
  if (fimo2$orientation[row] == "H-H") {
    fimo2$motif2[row] <- as.character(reverseComplement(DNAString(fimo2$motif2[row])))
  } else if (fimo2$orientation[row] == "T-T") {
    fimo2$motif1[row] <- as.character(reverseComplement(DNAString(fimo2$motif1[row])))
  }
}


#calculate # mismatches and categorize
fimo2$motif1mis <- as.vector(t(adist("TGAAACA", fimo2$motif1)))
fimo2$motif2mis <- as.vector(t(adist("TGAAACA", fimo2$motif2)))

#combine 0-1 and 1-0, etc. as one group
#not sure how to do all at once but easy if looping
fimo2$motifcat <- NA
for (row in 1:nrow(fimo2)) {
  if (!(is.na(fimo2$motif1mis[row]))) {
    if (fimo2$motif1mis[row] < fimo2$motif2mis[row]) {
      fimo2$motifcat[row] <- paste(fimo2$motif1mis[row], fimo2$motif2mis[row], sep="_")
    } else if (fimo2$motif1mis[row] >= fimo2$motif2mis[row]) {
      fimo2$motifcat[row] <- paste(fimo2$motif2mis[row], fimo2$motif1mis[row], sep="_")
    }
  }
}

fimo2$mismatchsum <- rowSums(fimo2[,c("motif1mis", "motif2mis")])

# View(ftable(fimo2$motif1 ~ fimo2$motif2))
# View(ftable(fimo2$motifcat))
# View(ftable(fimo2$mismatchsum))


#merge with peaks information to look at log2 FC, Q value, etc.
peaks2$shortname <- sapply(peaks2$newname, FUN = function(x){strsplit(x, "_")[[1]][1]})

#first have to fix fimo2 names by replacing first _ with a :
#luckily since I only want first occurrence, can use sub
fimo2$sequence_name2 <- sapply(fimo2$sequence_name, FUN = function(x){sub(x, pattern = "_", replacement = ":")})
fimo2$shortname <- sapply(fimo2$sequence_name2, FUN = function(x){strsplit(x, "_")[[1]][1]})

peaks3 <- merge(peaks2, fimo2, by="shortname", all=T)

#save final file
write.table(peaks3, file="analysis/1_1_Kar4dep_Ste12_AF.txt_win24/1_1_Kar4dep_Ste12_AF.txt_win24_1_1_Kar4dep_Ste12_AF.events.curated.txt_dupsmerged_241bp_nc.fasta_fimo_curated.txt", sep="\t", quote=F, row.names = F)








