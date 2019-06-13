#lodscore curves

setwd("C:/Users/Manuel/Desktop/AIL_B6xBFMI/RAWDATA")

#genotypes <- read.csv("genomatrix.clean.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
genotypes <- read.csv("genomatrix.clean.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
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

# Make sure that the ordering between phenotypes and genotypes matches !!!!!
# Also sort the markers by their natural chromosome ordering
genotypes <- genotypes[rownames(annotation), rownames(phenotypes)]

## Chr 3 bodyweight
# Dataset, containing columns named: Chr, Pos, marPvalue
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

# Day 63 with effect plot
par(mfrow=c(2,2))
plot(main = "Lod score curve Chr 3 (day 63)", c(min(as.numeric(chr3[, "Position"])), max(as.numeric(chr3[, "Position"]))), c(0,9), ylab = "Lodscore", xlab = "Position")
	points(x = as.numeric(chr3[,"Position"]), y = chr3[,"d63"] , type = "l", col="deepskyblue", lwd = 1)
	abline(h=4.5, col="red")

topmarkerID <- rownames(dataset[which.max(dataset[,"d63"]),])
topmarker <- t(genotypes[topmarkerID,])
groupsSize <- apply(genotypes,1,  table)[[topmarkerID]]
genopheno <- cbind(topmarker, phenotypes[,"d63"])
colnames(genopheno) <- c("Genotype", "d63")
bpt <- boxplot(as.numeric(as.character(genopheno[, "d63"]))  ~ genopheno[,"Genotype"], main = "d63 Topmarker", xlab = "Genotype", ylab = "weight (gr.)", col = c("lightblue", "lightgreen", "pink"),  notch = TRUE)
lines(1:3, bpt$stats[ 3, ], col="red", lwd=1)
	legend("bottomright", 
	legend = c(groupsSize[1], groupsSize[2], groupsSize[3]), 
       col = c("lightblue", "lightgreen", "pink"),
       pch = c(19, 19), cex = 0.8, 
       box.col = "darkgreen"
       )
	 
# Day 98 with effect plot
plot(main = "Lod score curve Chr 3 (day 98)", c(min(as.numeric(chr3[, "Position"])), max(as.numeric(chr3[, "Position"]))), c(0,9), ylab = "Lodscore", xlab = "Position")
	points(x = as.numeric(chr3[,"Position"]), y = chr3[,"d98"] , type = "l", col="deepskyblue", lwd = 1)
	abline(h=4.5, col="red")

topmarkerID <- rownames(dataset[which.max(dataset[,"d98"]),])
topmarker <- t(genotypes[topmarkerID,])
groupsSize <- apply(genotypes,1,  table)[[topmarkerID]]
genopheno <- cbind(topmarker, phenotypes[,"d98"])
colnames(genopheno) <- c("Genotype", "d98")
bpt <- boxplot(as.numeric(as.character(genopheno[, "d98"]))  ~ genopheno[,"Genotype"], main = "d98 Topmarker", xlab = "Genotype", ylab = "weight (gr.)", col = c("lightblue", "lightgreen", "pink"), notch = TRUE)
lines(1:3, bpt$stats[ 3, ], col="red", lwd=1)
	legend("bottomright", 
	legend = c(groupsSize[1], groupsSize[2], groupsSize[3]), 
       col = c("lightblue", "lightgreen", "pink"),
       pch = c(19, 19), cex = 0.8, 
       box.col = "darkgreen"
       )
	 
#Chr 8
#Dataset, containing columns named: Chr, Pos, marPvalue
dataset <- lodannotmatrix[, c("Chromosome", "Position", "Triglycerides/Proteins")]
chr8 <- dataset[which(dataset[,"Chromosome"] == 8),]

plot(main = "Lod score curve Chr 8", c(min(as.numeric(chr8[, "Position"])), max(as.numeric(chr8[, "Position"]))), c(0,7), ylab = "Lodscore", xlab = "Position")
	points(x = as.numeric(chr8[,"Position"]), y = chr8[,"Triglycerides/Proteins"] , type = "l", col="deepskyblue", lwd = 1)
	abline(h=4.5, col="red")
	
topmarkerID <- rownames(dataset[which.max(dataset[,"Triglycerides/Proteins"]),])
topmarker <- t(genotypes[topmarkerID,])
groupsSize <- apply(genotypes,1,  table)[[topmarkerID]]
genopheno <- cbind(topmarker, phenotypes[,"Triglycerides/Proteins"])
colnames(genopheno) <- c("Genotype", "Triglycerides/Proteins")
bpt <- boxplot(as.numeric(as.character(genopheno[, "Triglycerides/Proteins"]))  ~ genopheno[,"Genotype"], main = "Triglycerides/Proteins Topmarker", xlab = "Genotype", ylab = "Triglycerides/Proteins", col = c("lightblue", "lightgreen", "pink"), notch = TRUE)
lines(1:3, bpt$stats[ 3, ], col="red", lwd=1)
	legend("bottomright", 
	legend = c(groupsSize[1], groupsSize[2], groupsSize[3]), 
       col = c("lightblue", "lightgreen", "pink"),
       pch = c(19, 19), cex = 0.8, 
       box.col = "darkgreen"
       )