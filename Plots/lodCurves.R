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
dataset <- lodannotmatrix[, c("Chromosome", "Position", "d63", "d70", "d98", "d105", "d140")]
chr3 <- dataset[which(dataset[,"Chromosome"] == 3),]
signMarkers <- dataset[which(dataset[,"d98"] > 4.5),]

# lodcurve plots for each time point
plot(main = "Lod score curve Chr 3", c(min(as.numeric(chr3[, "Position"])), max(as.numeric(chr3[, "Position"]))), c(0,9), ylab = "-log10(pvalue)", xlab = "Position", las = 2, t = "n", xaxt = "n")
  points(x = as.numeric(chr3[,"Position"]), y = chr3[,"d70"] , type = "l", col="dodgerblue", lwd = 1)
  points(x = as.numeric(chr3[,"Position"]), y = chr3[,"d140"] , type = "l", col="blue", lwd = 1)
  points(x = as.numeric(chr3[,"Position"]), y = chr3[,"d63"] , type = "l", col="deepskyblue", lwd = 1)
  points(x = as.numeric(chr3[,"Position"]), y = chr3[,"d105"] , type = "l", col="purple", lwd = 1)
  points(x = as.numeric(chr3[,"Position"]), y = chr3[,"d98"] , type = "l", col="dodgerblue4", lwd = 1)
  abline(h=4.5, col="green")
  abline(h=4, col="orange")
  axis(1, at = c(0,25000000, 50000000, 75000000, 100000000, 125000000, 150000000), c("0 mb", "25 mb", "50 mb", "75 mb", "100 mb", "125 mb", "150 mb"))
  legend("topleft",
   legend = c("d63", "d70", "d98", "d105", "d140"),
    col = c("deepskyblue", "dodgerblue", "dodgerblue4", "blue", "purple"),
    pch = c(20,20,20),
    bty = "n",
    pt.cex = 2,
    cex = 1.2,
    text.col = "black")
	 
# Growth curve for each timepoint using the main topmarker
topmarkerID <- rownames(dataset[which.max(dataset[,"d63"]),])
topmarker <- t(genotypes[topmarkerID,])
groupsSize <- apply(genotypes,1,  table)[[topmarkerID]]
genopheno <- cbind(topmarker, phenotypes[,"d63"], phenotypes[,"d70"], phenotypes[,"d98"], phenotypes[,"d105"], phenotypes[,"d140"])
genopheno <- gsub("AA" ,"BFMI", genopheno)
genopheno <- gsub("AG" ,"HET", genopheno)
genopheno <- gsub("GG", "B6", genopheno)
colnames(genopheno) <- c("Genotype", 63, 70, 98, 105, 140)	 
timepoints <- as.numeric(as.character((colnames(genopheno[,-1]))))
plot(main="growth curve", c(60,140), c(0,70), ylab="weight (grams)", xlab="Time points", las = 2, t = "n", xaxt="n")
axis(1, at = c(63, 70, 98, 105, 140), c("day 60", "day 70", "day 98", "day 105", "day 140"))
meansBFMI <- c()
meansHET <- c()
meansB6 <- c()
for (x in timepoints){
 bptBFMI <- boxplot(at = x+1, as.numeric(as.character(genopheno[which(genopheno == "BFMI"), as.character(x)]))  ~ genopheno[which(genopheno == "BFMI"),"Genotype"], col = "lightgreen", axes = FALSE, add=TRUE) 
 meanBFMI <- bptBFMI$stats[3,] 
 meansBFMI <- c(meansBFMI, meanBFMI)
 bptHET <- boxplot(at = x, as.numeric(as.character(genopheno[which(genopheno == "HET"), as.character(x)]))  ~ genopheno[which(genopheno == "HET"),"Genotype"], col = "pink", axes = FALSE, add=TRUE) 
 meanHET <- bptHET$stats[3,]
 meansHET <- c(meansHET, meanHET)
 bptB6 <- boxplot(at = x-1, as.numeric(as.character(genopheno[which(genopheno == "B6"), as.character(x)]))  ~ genopheno[which(genopheno == "B6"),"Genotype"], col = "lightblue", axes = FALSE, add=TRUE) 
 meanB6 <- bptB6$stats[3,]
 meansB6 <- c(meansB6, meanB6)
 }
lines(c(64, 71, 99, 106, 141), meansBFMI, col="green", lwd=1)
lines(c(63, 70, 98, 105, 140), meansHET, col="pink", lwd=1)
lines(c(62, 69, 97, 104, 139), meansB6, col="blue", lwd=1) 
legend("topright", 
 legend = c("B6", "BFMI", "HET"), 
 col = c("lightblue", "lightgreen", "pink"),
 pch = 15, cex = 0.8, 
 box.col = "darkgreen")
 
 
# Day 63 with effect plot
par(mfcol=c(2,2))
plot(main = "Lod score curve Chr 3 (day 63)", c(min(as.numeric(chr3[, "Position"])), max(as.numeric(chr3[, "Position"]))), c(0,9), ylab = "Lodscore", xlab = "Position")
	points(x = as.numeric(chr3[,"Position"]), y = chr3[,"63"] , type = "l", col="deepskyblue", lwd = 1)
	abline(h=4.5, col="red")
	abline(h=4, col="green")
	
bpt <- boxplot(as.numeric(as.character(genopheno[, "63"]))  ~ genopheno[,"Genotype"], main = "d63 Topmarker", xlab = "Genotype", ylab = "weight (gr.)", col = c("lightblue", "lightgreen", "pink"),  notch = TRUE)
	   lines(1:3, bpt$stats[ 3, ], col="red", lwd=1)
	   legend("bottomright", 
	   legend = c(groupsSize[1], groupsSize[2], groupsSize[3]), 
       col = c("lightblue", "lightgreen", "pink"),
       pch = c(19, 19), cex = 0.8, 
       box.col = "darkgreen"
       )
	 
# Day 98 with effect plot
plot(main = "Lod score curve Chr 3 (day 98)", c(min(as.numeric(chr3[, "Position"])), max(as.numeric(chr3[, "Position"]))), c(0,9), ylab = "Lodscore", xlab = "Position")
	points(x = as.numeric(chr3[,"Position"]), y = chr3[,"98"] , type = "l", col="deepskyblue", lwd = 1)
	abline(h=4.5, col="red")
	abline(h=4, col="green")

topmarkerID <- rownames(dataset[which.max(dataset[,"98"]),])
topmarker <- t(genotypes[topmarkerID,])
groupsSize <- apply(genotypes,1,  table)[[topmarkerID]]
genopheno <- cbind(topmarker, phenotypes[,"98"])
genopheno <- gsub("AA" ,"BFMI", genopheno)
genopheno <- gsub("AG" ,"HET", genopheno)
genopheno <- gsub("GG", "B6", genopheno)
colnames(genopheno) <- c("Genotype", "98")

bpt <- boxplot(as.numeric(as.character(genopheno[, "98"]))  ~ genopheno[,"Genotype"], main = "d98 Topmarker", xlab = "Genotype", ylab = "weight (gr.)", col = c("lightblue", "lightgreen", "pink"), notch = TRUE)
lines(1:3, bpt$stats[ 3, ], col="red", lwd=1)
	legend("topright", 
	legend = c(groupsSize[1], groupsSize[2], groupsSize[3]), 
       col = c("lightblue", "lightgreen", "pink"),
       pch = c(19, 19), cex = 0.8, 
       box.col = "darkgreen"
       )
	   
# Chr 8
# Dataset, containing columns named: Chr, Pos, marPvalue
dataset <- lodannotmatrix[, c("Chromosome", "Position", "Triglycerides/Proteins")]
chr8 <- dataset[which(dataset[,"Chromosome"] == 8),]

plot(main = "Lod score curve Chr 8", c(min(as.numeric(chr8[, "Position"])), max(as.numeric(chr8[, "Position"]))), c(0,7), ylab = "Lodscore", xlab = "Position")
	points(x = as.numeric(chr8[,"Position"]), y = chr8[,"Triglycerides/Proteins"] , type = "l", col="deepskyblue", lwd = 1)
	abline(h=4.5, col="red")
	abline(h=4, col="green")
	
topmarkerID <- rownames(dataset[which.max(dataset[,"Triglycerides/Proteins"]),])
topmarker <- t(genotypes[topmarkerID,])
groupsSize <- apply(genotypes,1,  table)[[topmarkerID]]
genopheno <- cbind(topmarker, phenotypes[,"Triglycerides/Proteins"])
colnames(genopheno) <- c("Genotype", "Triglycerides/Proteins")
bpt <- boxplot(as.numeric(as.character(genopheno[, "Triglycerides/Proteins"]))  ~ genopheno[,"Genotype"], main = "Triglycerides/Proteins Topmarker", xlab = "Genotype", ylab = "Triglycerides/Proteins", col = c("lightblue", "lightgreen", "pink"), notch = TRUE)
lines(1:3, bpt$stats[ 3, ], col="red", lwd=1)
	legend("topright", 
	legend = c(groupsSize[1], groupsSize[2], groupsSize[3]), 
       col = c("lightblue", "lightgreen", "pink"),
       pch = c(19, 19), cex = 0.8, 
       box.col = "darkgreen"
       )
	   
# Chr 1
# Dataset, containing columns named: Chr, Pos, marPvalue
dataset <- lodannotmatrix[, c("Chromosome", "Position", "LiverWeight")]
chr1 <- dataset[which(dataset[,"Chromosome"] == 1),]

plot(main = "Lod score curve Chr 1", c(min(as.numeric(chr1[, "Position"])), max(as.numeric(chr1[, "Position"]))), c(0,7), ylab = "Lodscore", xlab = "Position")
	points(x = as.numeric(chr1[,"Position"]), y = chr1[,"LiverWeight"] , type = "l", col="deepskyblue", lwd = 1)
	abline(h=4.5, col="red")
	abline(h=4, col="green")

mmodel <- lm(phenotypes[,"Leber"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
res <- residuals(mmodel)
correctPhen <- residuals + intercept
	
topmarkerID <- rownames(chr1[which.max(chr1[,"LiverWeight"]),])
topmarker <- t(genotypes[topmarkerID,])
groupsSize <- apply(genotypes,1, table)[[topmarkerID]]
genopheno <- cbind(topmarker, phenotypes[,"Leber"])
colnames(genopheno) <- c("Genotype", "LiverWeight")
bpt <- boxplot(as.numeric(as.character(genopheno[,"LiverWeight"]))  ~ genopheno[,"Genotype"], main = "LiverWeight Topmarker", xlab = "Genotype", ylab = "LiverWeight", col = c("lightblue", "lightgreen", "pink"), notch = TRUE)
lines(1:2, bpt$stats[ 3, ], col="red", lwd=1)
	legend("topright", 
	legend = c(groupsSize[1], groupsSize[2]), 
       col = c("lightblue", "lightgreen"),
       pch = c(19, 19), cex = 0.8, 
       box.col = "darkgreen"
       )
