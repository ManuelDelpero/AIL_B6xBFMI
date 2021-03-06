setwd("C:/Users/Manuel/Desktop/AIL_B6xBFMI/RAWDATA")

rawgeno <- read.csv("RawGenotypes.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses=c(rep("character",4), rep("numeric",7)), na.strings=c("-", "NaN", ""))
markernames <- unique(rawgeno[,"SNP Name"])
samplenames <- gsub("666", "666 ", unique(rawgeno[,"Sample ID"]))

genotypes <- matrix(NA, length(markernames), length(samplenames), dimnames=list(markernames, samplenames))
for(x in 1:nrow(rawgeno)){
  if(!is.na(rawgeno[x, "GC Score"]) && rawgeno[x, "GC Score"] > 0.7) {
    sname <- gsub("666", "666 ", unique(rawgeno[x,"Sample ID"]))
    genotypes[rawgeno[x,"SNP Name"], sname] <- paste0(rawgeno[x,c("Allele1 - Forward", "Allele2 - Forward")], collapse="")
  }
}

# We could do some additional visual checks
# plot(rawgeno[which(rawgeno[,1] == "B10010006645"),"X"],rawgeno[which(rawgeno[,1] == "B10010006645"),"Y"])
# plot(rawgeno[which(rawgeno[,1] == "B10010021636"),"X"],rawgeno[which(rawgeno[,1] == "B10010021636"),"Y"])

write.table(genotypes, "genomatrix.raw.txt", sep = "\t", quote=FALSE)

# Dimensions we start with
dim(genotypes)

# QC on genotypes
# All genotypes are missing
idx <- which(apply(apply(genotypes, 1, is.na), 2, sum) == ncol(genotypes))
genotypes <- genotypes[-idx, ]
dim(genotypes)

# genotype is not segregating
nonSeg <- which(unlist(lapply(apply(genotypes,1,table), function(x){length(x) == 1})))
genotypes <- genotypes[-nonSeg, ]
dim(genotypes)

# At least 2 groups with 10 observations
good <- which(unlist(lapply(apply(genotypes,1,table), function(x){
  length(which(x > 10)) >= 2
})))
genotypes <- genotypes[good, ]
dim(genotypes)

# Prevent small groups, set the groups with < 10 individuals to NA
mtab <- apply(genotypes,1,table)
for(x in 1:nrow(genotypes)){
  if(any(mtab[[x]] < 6)){
    gt <- names(mtab[[x]])[which(mtab[[x]] < 6)]
    genotypes[x, genotypes[x, ] == gt] <- NA
  }
}
dim(genotypes)

# No duplicated markers, keep the first one we see
genotypes <- genotypes[-which(duplicated(genotypes)),]
dim(genotypes)

write.table(genotypes, "genomatrix.clean.txt", sep = "\t", quote=FALSE)

# Conver genotypes to numerical values to map using an additive model
numgeno <- matrix(NA, nrow(genotypes), ncol(genotypes), dimnames=list(rownames(genotypes), colnames(genotypes)))
for(x in 1:nrow(genotypes)){
  alleles <- na.omit(unique(unlist(strsplit(genotypes[x,], ""))))
  h1 <- paste0(alleles[1], alleles[1])
  het <- c(paste0(alleles[1], alleles[2]), paste0(alleles[2], alleles[1]))
  h2 <- paste0(alleles[2], alleles[2])
  numgeno[x, which(genotypes[x, ] == h1)] <- -1
  numgeno[x, which(genotypes[x, ] %in% het)] <- 0
  numgeno[x, which(genotypes[x, ] == h2)] <- 1
}
write.table(numgeno, "genomatrix.clean_numeric.txt", sep = "\t", quote=FALSE)
