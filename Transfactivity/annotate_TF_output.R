
library(reshape2)

##have to call this before changing the working directory for easy transferring between systems
sysNames <- read.delim(file="SGD_features.tab", header=TRUE, stringsAsFactors=FALSE)

fillMissingNames <- function(data, base, replace) {
  i <- sapply(data, is.factor) #replace factors with characters
  data[i] <- lapply(data[i], as.character)
  for (name in 1:length(data[,base])) {
    if (data[name,base]=="" | is.na(data[name,base])==TRUE) {
      data[name,base] <- data[name,replace]
    }
  }
  return(data)
}



###note that this function does not perfectly preserve row order. Do not assume it will.
swapSystematicForCommon <- function(data, systematic, names) { #data refers to data matrix, systematic refers to column with ORFs
  data2 <- merge(data, names, by.x=systematic, by.y=4, all.x=TRUE, sort=FALSE) #big merge
  startCol <- length(data[1,])
  data3 <- data2[,c(1:startCol, startCol+4)] #isolate just gene names
  data3 <- fillMissingNames(data3, length(data3), 1) #fill missing names
  data3 <- data3[,c(length(data3), 2:(length(data3)-1))] #Always puts gene names as first column now
  return(data3)
}


args2 <- commandArgs(trailingOnly=TRUE)

#args are working directory and input directory
args <- c("",
          "transfactivity_mvar_coeff.tsv")

if (length(args2)>0) {
  for (i in 1:length(args2)) {
    args[i] <- args2[i]
  }
}


setwd(args[1])
input_file <- args[2]
input_matrix <- read.delim(file=input_file, header=FALSE, stringsAsFactors = FALSE)

#transpose and remove extra info
input_matrix <- t(input_matrix)
input_matrix <- input_matrix[c(1,4:nrow(input_matrix)),]


#now convert p-values, if working on a pval matrix
if (args[2]=="transfactivity_mvar_pvals.tsv") {
  #convert to -log10
  if (ncol(input_matrix)>2) {
    input_matrix[2:nrow(input_matrix),2:ncol(input_matrix)] <- apply(input_matrix[2:nrow(input_matrix),2:ncol(input_matrix)],2,function(x) {-log10(as.numeric(x))})
  } else {
    input_matrix[2:nrow(input_matrix),2:ncol(input_matrix)] <- -log10(as.numeric(input_matrix[2:nrow(input_matrix),2:ncol(input_matrix)]))
  }
  #get sign of coefficient
  tmpFile <- read.delim(file="transfactivity_mvar_coeff.tsv", header=FALSE, stringsAsFactors = FALSE)
  tmpFile <- t(tmpFile)
  tmpFile <- tmpFile[c(1,4:nrow(tmpFile)),]
  #now get the sign and multiply by the p-value matrix
  if (ncol(input_matrix)>2) {
    input_matrix[2:nrow(input_matrix), 2:ncol(input_matrix)] <- as.numeric(input_matrix[2:nrow(input_matrix), 2:ncol(input_matrix)])*apply(tmpFile[2:nrow(tmpFile), 2:ncol(tmpFile)], 2, function(x) {sign(as.numeric(x))})
  } else {
    input_matrix[2:nrow(input_matrix), 2:ncol(input_matrix)] <- as.numeric(input_matrix[2:nrow(input_matrix), 2:ncol(input_matrix)])*sign(as.numeric(tmpFile[2:nrow(tmpFile), 2:ncol(tmpFile)]))
  }
}

#now map PSAM ID back to TF ID
PSAMs <- read.delim(file="design_matrix.tsv", header=FALSE, stringsAsFactors = FALSE, nrows=nrow(input_matrix)-1, skip=1)

#now parse out the PSAM strings, split on "/" and take the final value
PSAMlist <- apply(PSAMs, 1, function(x) {strsplit(x, split="/")[[1]][length(strsplit(x, split="/")[[1]])]})

#now remove the ".xml" string
PSAMlist <- sub(pattern=".xml", x=PSAMlist, replacement="")

#now update the values in the original matrix. This assumes that everything is still in the same order.... should almost certainly be true
input_matrix[2:nrow(input_matrix),1] <- PSAMlist

#convert top row to column names
colnames(input_matrix) <- input_matrix[1,]
input_matrix <- input_matrix[2:nrow(input_matrix),]


######swap systematic names for common section. Needlessly complicated but it works
geneList <- strsplit(input_matrix[,1], "_")

#build data frame
geneList2 <- vector()
for (gene in 1:length(geneList)) {
  geneList2 <- c(geneList2, geneList[[gene]][1])
}

#put it in a data frame because that's what swap function expects
genedf <- data.frame(geneList2, input_matrix)

genedf <- swapSystematicForCommon(genedf, 1, sysNames)

#now append the motif identifier
for (row in 1:nrow(genedf)) {
  append <- sub(x=strsplit(genedf[row,2], "_")[[1]][2], pattern=".pfm", replacement="")
  genedf[row,1] <- paste(genedf[row,1], append, sep="_")
}

#now remove the systematic column
genedf <- genedf[,c(1,3:ncol(genedf))]



#now save the result
write.table(genedf, file=paste0(input_file,"_annotated.txt",sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)





