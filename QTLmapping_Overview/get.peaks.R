# Call peaks in the lodmatrix
#
# copyright (c) 2015-2020 - Brockmann group - HU Berlin, Danny Arends & Manuel DelPero
# last modified May, 2019
# first written May, 2019

setwd("C:/Users/Manuel/Desktop/AIL_B6xBFMI/RAWDATA")
source("C:/Github/AIL_B6xBFMI/QTLmapping_Overview/peak.detect.R")

genotypes <- read.table("OrderedGenotypes.txt", sep = "\t", check.names=FALSE)
markerannot <- read.csv("SNP_Map.txt", header=TRUE, sep="\t", row.names=2, check.names=FALSE)
markerannot <- markerannot[rownames(genotypes),]

phenotypes <- read.csv("allPhenotypes.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
mprofiles <- read.table("lodmatrix.adj.txt", sep="\t", row.names=1, header=TRUE)

results <- NULL
for(x in 1:ncol(mprofiles)){
  peaks <- peak.detect(mprofiles[,x], markerannot, loddrop = 1.5)
  if(!is.null(peaks)){
    for(p in 1:nrow(peaks)){
      leftPos <- markerannot[peaks[p,1], "Position"]
      topPos <- markerannot[peaks[p,2], "Position"]
      rightPos <- markerannot[peaks[p,3], "Position"]
      topChr <- as.character(markerannot[peaks[p,2], "Chromosome"])
      results <- rbind(results, c(colnames(mprofiles)[x], topChr, leftPos, topPos, rightPos, round(mprofiles[peaks[p,2],x],2), rownames(markerannot)[peaks[p,]]))
    }
  }
}

colnames(results) <- c("Phenotype", "Chr", "StartPos", "TopPos", "StopPos", "LOD", "flankLeft", "TopMarker", "FlankRight")
results <- data.frame(results)

write.table(results, "QTL_regions.txt", row.names=FALSE, quote=FALSE, sep='\t')
