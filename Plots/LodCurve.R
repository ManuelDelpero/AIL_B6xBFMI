#lodscore curves

setwd("C:/Users/Manuel/Desktop/AIL_B6xBFMI/RAWDATA")

#genotypes <- read.csv("genomatrix.clean.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
genotypes <- read.csv("genomatrix.clean_numeric.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
phenotypes <- read.csv("allPhenotypes.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
markerannot <- read.csv("SNP_Map.txt", header=TRUE, sep="\t", row.names=2, check.names=FALSE)
lodmatrix <- read.csv("lodmatrix.adj.txt", header=TRUE, sep="\t", check.names=FALSE)
markerannot <- markerannot[,-1]
markerannot <- markerannot[rownames(genotypes),]
markerannot <- markerannot[sort(markerannot[,"Position"], index.return=TRUE)$ix,]

chromosomes <- c(1:19, "X", "Y")

annotation <- c()
for(chr in chromosomes){
  annotation <- rbind(annotation, markerannot[markerannot[,"Chromosome"] == chr,])
}

lodannotmatrix <- cbind(annotation[rownames(lodmatrix), ], lodmatrix)




#Chr 3 bodyweight
#Dataset, containing columns named: Chr, Pos, marPvalue
dataset <- lodannotmatrix[, c("Chromosome", "Position", "d63", "d70", "d98", "d140")]
chr3 <- dataset[which(dataset[,"Chromosome"] == 3),]

plot(main = "Lod score curve Chr 3", c(min(as.numeric(chr3[, "Position"])), max(as.numeric(chr3[, "Position"]))), c(0,9), ylab = "Lodscore", xlab = "Position")
	points(x = as.numeric(chr3[,"Position"]), y = chr3[,"d70"] , type = "l", col="dodgerblue", lwd = 1)
	points(x = as.numeric(chr3[,"Position"]), y = chr3[,"d140"] , type = "l", col="blue", lwd = 1)
	points(x = as.numeric(chr3[,"Position"]), y = chr3[,"d63"] , type = "l", col="deepskyblue", lwd = 1)
	points(x = as.numeric(chr3[,"Position"]), y = chr3[,"d98"] , type = "l", col="dodgerblue4", lwd = 1)
	abline(h=4.5, col="red")
  legend("topleft",
  legend = c("d63", "d70", "d98", "d140"),
     col = c("deepskyblue", "dodgerblue", "dodgerblue4", "blue"),
     pch = c(20,20,20),
     bty = "n",
     pt.cex = 2,
     cex = 1.2,
     text.col = "black")



#Chr 8
#Dataset, containing columns named: Chr, Pos, marPvalue
dataset <- lodannotmatrix[, c("Chromosome", "Position", "Triglycerides/Proteins")]
chr3 <- dataset[which(dataset[,"Chromosome"] == 8),]

plot(main = "Lod score curve Chr 8", c(min(as.numeric(chr3[, "Position"])), max(as.numeric(chr3[, "Position"]))), c(0,7), ylab = "Lodscore", xlab = "Position")
	points(x = as.numeric(chr3[,"Position"]), y = chr3[,"Triglycerides/Proteins"] , type = "l", col="orange", lwd = 3)
	abline(h=4.5, col="red")
