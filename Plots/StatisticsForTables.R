# Basic statistics for manuscript tables
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

## Figure out means for each genotypes

# Triglycerides
data <- lodannotmatrix[, c("Chromosome", "Position", "Triglycerides/Proteins")]
chr8 <- data[which(Trig[,"Chromosome"] == 8),]
topmarkerID <- rownames(data[which.max(data[,"Triglycerides/Proteins"]),])
topmarker <- t(genotypes[topmarkerID,])
groupsSize <- apply(genotypes,1,  table)[[topmarkerID]]
genopheno <- cbind(topmarker, phenotypes[,"Triglycerides/Proteins"])
colnames(genopheno) <- c("Genotype", "Triglycerides/Proteins")
meanBFMI <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "AA"), 2]), na.rm = TRUE)
meanB6 <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "GG"), 2]), na.rm = TRUE)
meanHET <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "AG"), 2]), na.rm = TRUE)
meanBFMI
meanB6
meanHET

# Bodyweight chr 3
data <- lodannotmatrix[, c("Chromosome", "Position", "d140")]
chr3 <- data[which(data[,"Chromosome"] == 3),]
topmarkerID <- rownames(data[which.max(data[,3]),])
topmarker <- t(genotypes[topmarkerID,])
groupsSize <- apply(genotypes,1,  table)[[topmarkerID]]
genopheno <- cbind(topmarker, phenotypes[,"d140"])
colnames(genopheno) <- c("Genotype", "d63")
meanBFMI <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "AA"), 2]), na.rm = TRUE)
meanB6 <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "GG"), 2]), na.rm = TRUE)
meanHET <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "AG"), 2]), na.rm = TRUE)
meanBFMI
meanB6
meanHET

# Bodyweight chr 6
data <- lodannotmatrix[, c("Chromosome", "Position", "d112")]
data <- data[which(data[,"Chromosome"] == 6),]
topmarkerID <- rownames(data[which.max(data[,3]),])
topmarker <- t(genotypes[topmarkerID,])
groupsSize <- apply(genotypes,1,  table)[[topmarkerID]]
genopheno <- cbind(topmarker, phenotypes[,"d112"])
colnames(genopheno) <- c("Genotype", "d63")
meanBFMI <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "TT"), 2]), na.rm = TRUE)
meanB6 <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "GG"), 2]), na.rm = TRUE)
meanHET <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "TG"), 2]), na.rm = TRUE)
meanBFMI
meanB6
meanHET

# Liver weight
data <- lodannotmatrix[, c("Chromosome", "Position", "Leber")]
chr8 <- data[which(data[,"Chromosome"] == 1),]
topmarkerID <- rownames(data[which.max(data[,"Leber"]),])
topmarker <- t(genotypes[topmarkerID,])
groupsSize <- apply(genotypes,1,  table)[[topmarkerID]]
genopheno <- cbind(topmarker, phenotypes[,"Leber"])
colnames(genopheno) <- c("Genotype", "Leber")
meanBFMI <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "AA"), 2]), na.rm = TRUE)
meanB6 <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "GG"), 2]), na.rm = TRUE)
meanHET <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "AG"), 2]), na.rm = TRUE)
meanBFMI
meanB6
meanHET

# Sub fat weight
data <- lodannotmatrix[, c("Chromosome", "Position", "WATsc")]
chr8 <- data[which(data[,"Chromosome"] == 1),]
topmarkerID <- rownames(data[which.max(data[,"WATsc"]),])
topmarker <- t(genotypes[topmarkerID,])
groupsSize <- apply(genotypes,1,  table)[[topmarkerID]]
genopheno <- cbind(topmarker, phenotypes[,"WATsc"])
colnames(genopheno) <- c("Genotype", "WATsc")
meanBFMI <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "AA"), 2]), na.rm = TRUE)
meanB6 <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "GG"), 2]), na.rm = TRUE)
meanHET <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "AG"), 2]), na.rm = TRUE)
meanBFMI
meanB6
meanHET

# Glucose
data <- lodannotmatrix[, c("Chromosome", "Position", "ITT_15")]
chr8 <- data[which(data[,"Chromosome"] == 1),]
topmarkerID <- rownames(data[which.max(data[,"ITT_15"]),])
topmarker <- t(genotypes[topmarkerID,])
groupsSize <- apply(genotypes,1,  table)[[topmarkerID]]
genopheno <- cbind(topmarker, phenotypes[,"ITT_15"])
colnames(genopheno) <- c("Genotype", "ITT_15")
meanBFMI <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "AA"), 2]), na.rm = TRUE)
meanB6 <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "GG"), 2]), na.rm = TRUE)
meanHET <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "AG"), 2]), na.rm = TRUE)
meanBFMI
meanB6
meanHET

# BMI
data <- lodannotmatrix[, c("Chromosome", "Position", "BMI")]
chr8 <- data[which(data[,"Chromosome"] == 1),]
topmarkerID <- rownames(data[which.max(data[,"BMI"]),])
topmarker <- t(genotypes[topmarkerID,])
groupsSize <- apply(genotypes,1,  table)[[topmarkerID]]
genopheno <- cbind(topmarker, phenotypes[,"BMI"])
colnames(genopheno) <- c("Genotype", "BMI")
meanBFMI <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "TT"), 2]), na.rm = TRUE)
meanB6 <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "CC"), 2]), na.rm = TRUE)
meanHET <- mean(as.numeric(genopheno[which(genopheno[, "Genotype"] == "TC"), 2]), na.rm = TRUE)
meanBFMI
meanB6
meanHET

