# Plots ideas for the manuscript
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin, Manuel Delpero
# last modified September, 2019
# first written August, 2019

setwd("C:/Users/Manuel/Desktop/AIL_B6xBFMI/RAWDATA")

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
dataset <- lodannotmatrix[, c("Chromosome", "Position", "d28", "d49", "d35", "d42", "d49", "d56", "d63", "d70", "d77", "d84", "d91", "d98", "d105", "d112", "d119", "d126", "d133", "d140")]
chr3 <- dataset[which(dataset[,"Chromosome"] == 3),]
signMarkers <- dataset[which(dataset[,"d98"] > 4.5),]

par(cex.lab=1.5, cex.main = 1.8, cex.axis = 1.6)
mat <- matrix(c(1,1,2,3), 2, 2, byrow = TRUE)
layout(mat, widths = rep.int(3, ncol(mat)))
# lodcurve plots for each time point
plot(main = "QTL profile bodyweight [Chr 3]", c(min(as.numeric(chr3[, "Position"])), max(as.numeric(chr3[, "Position"]))), c(0,9), ylab = "-log10 [pvalue]", xlab = "Position [mb]", las = 2, t = "n", xaxt = "n")
  points(x = as.numeric(chr3[,"Position"]), y = chr3[,"d70"] , type = "l", col="dodgerblue", lwd = 1)
  points(x = as.numeric(chr3[,"Position"]), y = chr3[,"d140"] , type = "l", col="blue", lwd = 1)
  points(x = as.numeric(chr3[,"Position"]), y = chr3[,"d63"] , type = "l", col="deepskyblue", lwd = 1)
  points(x = as.numeric(chr3[,"Position"]), y = chr3[,"d105"] , type = "l", col="purple", lwd = 1)
  points(x = as.numeric(chr3[,"Position"]), y = chr3[,"d98"] , type = "l", col="dodgerblue4", lwd = 1)
  abline(h=4.5, col="green")
  abline(h=4, col="orange")
  axis(1, at = c(0,25000000, 50000000, 75000000, 100000000, 125000000, 150000000), c("0", "25", "50", "75", "100", "125", "150"))
  legend("topright", bg="gray",
  legend = c("Week 9", "Week 10", "Week 14", "Week 15","Week 20"),
    col = c("deepskyblue", "dodgerblue", "dodgerblue4", "purple", "blue"),
    pch = 15,
    pt.cex = 1.7,
    pt.bg = "lightsteelblue1",
    cex = 1,
    text.col = "black")
 
# Growth curve for each timepoint using the main topmarker
mmodelD28 <- lm(phenotypes[,"d28"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
d28.adj <- resid(mmodelD28) + coef(mmodelD28)["(Intercept)"]
mmodelD35 <- lm(phenotypes[,"d35"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
d35.adj <- resid(mmodelD35) + coef(mmodelD35)["(Intercept)"]
mmodelD42 <- lm(phenotypes[,"d42"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
d42.adj <- resid(mmodelD42) + coef(mmodelD42)["(Intercept)"]
mmodelD49 <- lm(phenotypes[,"d49"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
d49.adj <- resid(mmodelD49) + coef(mmodelD49)["(Intercept)"]
mmodelD56 <- lm(phenotypes[,"d56"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
d56.adj <- resid(mmodelD56) + coef(mmodelD56)["(Intercept)"]
mmodelD63 <- lm(phenotypes[,"d63"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
d63.adj <- resid(mmodelD63) + coef(mmodelD63)["(Intercept)"]
mmodelD70 <- lm(phenotypes[,"d70"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
d70.adj <- resid(mmodelD70) + coef(mmodelD70)["(Intercept)"]
mmodelD77 <- lm(phenotypes[,"d77"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
d77.adj <- resid(mmodelD77) + coef(mmodelD77)["(Intercept)"]
mmodelD84 <- lm(phenotypes[,"d84"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
d84.adj <- resid(mmodelD84) + coef(mmodelD84)["(Intercept)"]
mmodelD91 <- lm(phenotypes[,"d91"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
d91.adj <- resid(mmodelD91) + coef(mmodelD91)["(Intercept)"]
mmodelD98 <- lm(phenotypes[,"d98"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
d98.adj <- resid(mmodelD98) + coef(mmodelD98)["(Intercept)"]
mmodelD105 <- lm(phenotypes[,"d105"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
d105.adj <- resid(mmodelD105) + coef(mmodelD105)["(Intercept)"]
mmodelD112 <- lm(phenotypes[,"d112"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
d112.adj <- resid(mmodelD112) + coef(mmodelD112)["(Intercept)"]
mmodelD119 <- lm(phenotypes[,"d119"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
d119.adj <- resid(mmodelD119) + coef(mmodelD119)["(Intercept)"]
mmodelD126 <- lm(phenotypes[,"d126"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
d126.adj <- resid(mmodelD126) + coef(mmodelD126)["(Intercept)"]
mmodelD133 <- lm(phenotypes[,"d133"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
d133.adj <- resid(mmodelD133) + coef(mmodelD133)["(Intercept)"]
mmodelD140 <- lm(phenotypes[,"d140"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
d140.adj <- resid(mmodelD140) + coef(mmodelD140)["(Intercept)"]
topmarkerID <- rownames(dataset[which.max(dataset[,"d63"]),])
topmarker <- t(genotypes[topmarkerID,])
groupsSize <- apply(genotypes,1,  table)[[topmarkerID]]
genopheno <- cbind(topmarker, d28.adj, d35.adj, d42.adj, d49.adj, d56.adj, d63.adj, d70.adj, d77.adj, d84.adj, d91.adj, d98.adj, d105.adj, d112.adj, d119.adj, d126.adj, d133.adj, d140.adj)
genopheno <- gsub("AA" ,"BFMI", genopheno)
genopheno <- gsub("AG" ,"HET", genopheno)
genopheno <- gsub("GG", "B6", genopheno)
colnames(genopheno) <- c("Genotype", 28, 35, 42, 49, 55, 63, 70, 77, 84, 91, 98, 105, 112, 119, 126, 133, 140)	 
timepoints <- as.numeric(as.character((colnames(genopheno[,-1]))))
timepoints <- c(28, 42, 55, 70, 84, 98, 112, 126, 140)  
 
plot(main="Body mass [genotype gUNC5036315]", c(25,140), c(0,55), ylab="Body mass [g.]", xlab="Age [weeks]", yaxs = "i", las = 2, t = "n", xaxt="n")
  #rect(50, 0, 102, 69.89, border = NA, col = "lightskyblue1")
  #rect(102, 0, 145, 69.89, border = NA, col = "lightskyblue")
  axis(1, at = c(28, 42, 56, 70, 84, 98, 112, 126, 140), c("4", "6", "8", "10", "12", "14", "16", "18", "20"), lwd = 1, cex.axis=1.6)
  meansBFMI <- c()
  meansHET <- c()
  meansB6 <- c()
  for (x in timepoints){
   bptBFMI <- boxplot(at = x+1, as.numeric(as.character(genopheno[which(genopheno == "BFMI"), as.character(x)]))  ~ genopheno[which(genopheno == "BFMI"),"Genotype"], width=6, col = "lightskyblue1", axes = FALSE, add=TRUE, notch= TRUE, outcex=0.5, boxwex = 1.5) 
   meanBFMI <- bptBFMI$stats[3,] 
   meansBFMI <- c(meansBFMI, meanBFMI)
   bptHET <- boxplot(at = x, as.numeric(as.character(genopheno[which(genopheno == "HET"), as.character(x)]))  ~ genopheno[which(genopheno == "HET"),"Genotype"], width=6, col = "cyan3", axes = FALSE, add=TRUE, notch= TRUE, outcex=0.5, boxwex = 1.5) 
   meanHET <- bptHET$stats[3,]
   meansHET <- c(meansHET, meanHET)
   bptB6 <- boxplot(at = x-1, as.numeric(as.character(genopheno[which(genopheno == "B6"), as.character(x)]))  ~ genopheno[which(genopheno == "B6"),"Genotype"], width=6, col = "dodgerblue4", axes = FALSE, add=TRUE, notch= TRUE, outcex=0.5, boxwex = 1.5) 
   meanB6 <- bptB6$stats[3,]
   meansB6 <- c(meansB6, meanB6)
  }
  lines(c(28, 42, 56, 70, 84, 98, 112, 126, 140), meansBFMI, col="lightskyblue1", lwd=1)
  lines(c(28, 42, 56, 70, 84, 98, 112, 126, 140), meansHET, col="cyan3", lwd=1)
  lines(c(28, 42, 56, 70, 84, 98, 112, 126, 140), meansB6, col="dodgerblue4", lwd=1) 
  legend("topleft",  
  legend = c("BFMI", "HET", "B6"), bg="gray", 
   col = c("lightskyblue1", "cyan3", "dodgerblue4"),
   pch = 15, pt.cex = 1.7, cex = 1, 
   box.col = "black")

# Day 105 effect plot for the second peak
topmarker <- t(genotypes["SBR032134332",])
groupsSize <- apply(genotypes,1,  table)[["SBR032134332"]]
genopheno <- cbind(topmarker, phenotypes[,"d105"])
genopheno <- gsub("AA" ,"BFMI", genopheno)
genopheno <- gsub("AG" ,"HET", genopheno)
genopheno <- gsub("GG", "B6", genopheno)
colnames(genopheno) <- c("Genotype", "d105")
bpt <- boxplot(as.numeric(as.character(genopheno[which(genopheno[,1] == "TT"),])), as.numeric(as.character(genopheno[which(genopheno[,1] == "TC"),])), as.numeric(as.character(genopheno[which(genopheno[,1] == "CC"),])), width = c(0.2, 0.2, 0.2), main = "Week 15 [genotype SBR032134332]", xlab = "Genotype", ylab = "Weight [g.]", col = c("dodgerblue4", "cyan3", "lightskyblue1"), las =2, xaxt = "n", ylim = c(0, 70))
 axis(1, at = 1:3 , c("CC", "TC", "TT"))
 lines(1:3, bpt$stats[3, ], col="red", lwd=1)
 legend("topleft", bg="gray", 
 legend = c("BFMI", "HET", "B6"), 
  col = c("lightskyblue1", "cyan3", "dodgerblue4"),
  pch = 15, pt.cex = 1.7, cex = 1, 
  box.col = "black")
 

# Day 63 with effect plot
#par(mfcol=c(2,2))
#plot(main = "Lod score curve Chr 3 (day 63)", c(min(as.numeric(chr3[, "Position"])), max(as.numeric(chr3[, "Position"]))), c(0,9), ylab = "Lodscore", xlab = "Position")
	#points(x = as.numeric(chr3[,"Position"]), y = chr3[,"63"] , type = "l", col="deepskyblue", lwd = 1)
	#abline(h=4.5, col="red")
	#abline(h=4, col="green")
	
#bpt <- boxplot(as.numeric(as.character(genopheno[, "63"]))  ~ genopheno[,"Genotype"], main = "d63 Topmarker", xlab = "Genotype", ylab = "weight (gr.)", col = c("lightblue", "lightgreen", "pink"),  notch = TRUE)
	   #lines(1:3, bpt$stats[ 3, ], col="red", lwd=1)
	   #legend("bottomright", 
	   #legend = c(groupsSize[1], groupsSize[2], groupsSize[3]), 
       #col = c("lightblue", "lightgreen", "pink"),
       #pch = c(19, 19), cex = 0.8, 
       #box.col = "darkgreen"
       #)
	 
# Day 98 with effect plot
#plot(main = "Lod score curve Chr 3 (day 98)", c(min(as.numeric(chr3[, "Position"])), max(as.numeric(chr3[, "Position"]))), c(0,9), ylab = "Lodscore", xlab = "Position")
	#points(x = as.numeric(chr3[,"Position"]), y = chr3[,"98"] , type = "l", col="deepskyblue", lwd = 1)
	#abline(h=4.5, col="red")
	#abline(h=4, col="green")

#topmarkerID <- rownames(dataset[which.max(dataset[,"98"]),])
#topmarker <- t(genotypes[topmarkerID,])
#groupsSize <- apply(genotypes,1,  table)[[topmarkerID]]
#genopheno <- cbind(topmarker, phenotypes[,"98"])
#genopheno <- gsub("AA" ,"BFMI", genopheno)
#genopheno <- gsub("AG" ,"HET", genopheno)
#genopheno <- gsub("GG", "B6", genopheno)
#colnames(genopheno) <- c("Genotype", "98")

#bpt <- boxplot(as.numeric(as.character(genopheno[, "98"]))  ~ genopheno[,"Genotype"], main = "d98 Topmarker", xlab = "Genotype", ylab = "weight (gr.)", col = c("lightblue", "lightgreen", "pink"), notch = TRUE)
#lines(1:3, bpt$stats[ 3, ], col="red", lwd=1)
	#legend("topright", 
	#legend = c(groupsSize[1], groupsSize[2], groupsSize[3]), 
       #col = c("lightblue", "lightgreen", "pink"),
       #pch = c(19, 19), cex = 0.8, 
       #box.col = "darkgreen"
       #)
# Chr 8
# Dataset, containing columns named: Chr, Pos, marPvalue

par(cex.lab=1.5, cex.main = 1.8, cex.axis = 1.4)
par(mfrow = c(2,2))
dataset <- lodannotmatrix[, c("Chromosome", "Position", "Triglycerides/Proteins")]
chr8 <- dataset[which(dataset[,"Chromosome"] == 8),]


plot(main = "QTL profile liver triglycerides [Chr 8]", c(min(as.numeric(chr8[, "Position"])), max(as.numeric(chr8[, "Position"]))), c(0,7), ylab = "-log10[pvalue]", xlab = "Position [mb]", las = 2, t = "n", xaxt = "n")
  points(x = as.numeric(chr8[,"Position"]), y = chr8[,"Triglycerides/Proteins"] , type = "l", col="deepskyblue", lwd = 1)
  abline(h=4.5, col="green")
  abline(h=4, col="orange")
  axis(1, at = c(0,25000000, 50000000, 75000000, 100000000, 125000000, 150000000), c("0", "25", "50", "75", "100", "125", "150"))

topmarkerID <- rownames(dataset[which.max(dataset[,"Triglycerides/Proteins"]),])
topmarker <- t(genotypes[topmarkerID,])
groupsSize <- apply(genotypes,1,  table)[[topmarkerID]]
genopheno <- cbind(topmarker, phenotypes[,"Triglycerides/Proteins"])
colnames(genopheno) <- c("Genotype", "Triglycerides/Proteins")
bpt <- boxplot(as.numeric(as.character(genopheno[which(genopheno[,1] == "GG"),])), as.numeric(as.character(genopheno[which(genopheno[,1] == "AG"),])), as.numeric(as.character(genopheno[which(genopheno[,1] == "AA"),])), main = "Liver triglycerides/Proteins[genotype S1H083826428]", xlab = "Genotype", ylab = "Triglycerides/Proteins", col = c("lightskyblue1", "cyan3", "dodgerblue4"), las =2, xaxt = "n")
  axis(1, at = 1:3 , c("AA", "AG", "GG"))
  lines(1:3, bpt$stats[ 3, ], col="red", lwd=1)
  legend("topleft", bg="gray", 
  legend = c("BFMI", "HET", "B6"), 
    col = c("lightskyblue1", "cyan3", "dodgerblue4"),
    pch = 15, pt.cex = 1.7, cex = 1, 
    box.col = "darkgreen")
 

# Chr 1
# Dataset, containing columns named: Chr, Pos, marPvalue
dataset <- lodannotmatrix[, c("Chromosome", "Position", "LiverWeight")]
chr1 <- dataset[which(dataset[,"Chromosome"] == 1),]

plot(main = "QTL profile liver weight [Chr 1]", c(min(as.numeric(chr1[, "Position"])), max(as.numeric(chr1[, "Position"]))), c(0,7), ylab = "-log10[pvalue]", xlab = "Position [mb]", las = 2, t = "n", xaxt = "n")
  points(x = as.numeric(chr1[,"Position"]), y = chr1[,"LiverWeight"] , type = "l", col="deepskyblue", lwd = 1)
  abline(h=4.5, col="green")
  abline(h=4, col="orange")
  axis(1, at = c(0,25000000, 50000000, 75000000, 100000000, 125000000, 150000000, 175000000, 200000000), c("0", "25", "50", "75", "100", "125", "150", "175", "200"))

mmodelLiver <- lm(phenotypes[,"Leber"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
liver.adj <- resid(mmodelLiver) + coef(mmodelLiver)["(Intercept)"]
	
topmarkerID <- rownames(chr1[which.max(chr1[,"LiverWeight"]),])
topmarker <- t(genotypes[topmarkerID,])
groupsSize <- apply(genotypes,1, table)[[topmarkerID]]
genopheno <- cbind(topmarker, liver.adj)
colnames(genopheno) <- c("Genotype", "LiverWeight")
bpt <- boxplot(as.numeric(as.character(genopheno[which(genopheno[,1] == "AG"),])), as.numeric(as.character(genopheno[which(genopheno[,1] == "GG"),])), main = "Liver Weight [genotype gUNC2036998]", xlab = "Genotype", ylab = "LiverWeight", col = c("lightskyblue1", "cyan3", "dodgerblue4"), las = 2, xaxt = "n")
  axis(1, at = 1:2 , c("AG", "GG"))
  lines(1:2, bpt$stats[ 3, ], col="red", lwd=1)
  legend("topleft", bg="gray",
  legend = c("HET", "BFMI"), 
    col = c("lightskyblue1", "cyan3"),
    pch = 15, pt.cex = 1.7, cex = 1, 
    box.col = "darkgreen")