setwd("C:\\Users\\Manuel\\Desktop\\miniMUGA")

                                                                       ## analysis on phenotypes ##


# read files
weight <- read.csv("weight.txt", header=TRUE, check.names=FALSE, sep="\t")
weight_test <- read.csv("weight_test.txt", header=TRUE, check.names=FALSE, sep="\t")
id <- read.csv("id_animals_Sex.txt", header=TRUE, check.names=FALSE, sep="\t")
id <- id[,-c(2,4)]
rownames(weight) <- id[,1]
weight_test <- weight_test[,-3]
weight <- cbind(weight, weight_test)
liv <- read.csv("Liv.txt", header=TRUE, check.names=FALSE, sep="\t")
rowMuscleFat <- read.csv("Musclefat.txt", header=TRUE, check.names=FALSE, sep="\t")
muscleFat <- rowMuscleFat[, c("Maus", "Fettanteil [%]")]
muscleFat <- muscleFat[-which(duplicated(muscleFat[,"Maus"])),]
rownames(muscleFat) <- gsub(" ", "", muscleFat[,1], fixed = TRUE)



#QC by weight
plot(main="weight growth", c(min(ncol(weight)), 140), c(0, max(weight, na.rm=TRUE)), t="n", ylim=c(0, max(weight, na.rm=TRUE)), yaxs="i", xlab="time(days)", ylab="weight(gr.)")
for (x in 1:nrow(weight))
 points(colnames(weight), weight[x,], t ="l", col = "green")

#weight lost oralGTT
plot(main="weight growth", c(100, 140), c(0, max(weight, na.rm=TRUE)), t="n", ylim=c(0, max(weight, na.rm=TRUE)), yaxs="i", xlab="time(days)", ylab="weight(gr.)")
for (x in 1:nrow(weight))
 points(colnames(weight), weight[x,], t ="l", col = "green")

# oralGTT
oralGTT <- read.csv("oralGTT.txt", header = TRUE, check.names = FALSE, sep = "\t")
time <- as.numeric(colnames(oralGTT))

plot(main =" oral GTT", c(min(time), max(time)), c(0, 600),  t = 'n', xlab="Time (min)", ylab="Blood Glucose (mg/dl)")
 for (row in 1:nrow(oralGTT))
   points(time, oralGTT[row,], t="l", col="blue")

bptOR <- boxplot(oralGTT, main="Oral Glucose Tolerance Test", ylab="Blood Glucose (mg/dl)", xlab="Time(min)", col=c("yellow"), notch=TRUE)
 lines(1:5, bptOR$stats[ 3, ], col="blue", lwd=2)

# insulinTT
insulinTT <- read.csv("insulinTT.txt", header=TRUE, check.names=FALSE, sep="\t")
time <- as.numeric(colnames(insulinTT))

plot(main =" insulinTT", c(min(time), max(time)), c(0, 600),  t = 'n', xlab="Time (min)", ylab="Blood Glucose (mg/dl)")
 for (y in 1:nrow(insulinTT))
   points(time, insulinTT[y,], t="l", col="blue")

bptIN <- boxplot(insulinTT, main="Insulin Tolerance Test", ylab="Blood Glucose (mg/dl)", xlab="Time(min)", col=c("yellow"), notch=TRUE)
 lines(1:4, bptIN$stats[ 3, ], col="blue", lwd=2)

#Tissues weight
tissuesWEIGHT <- read.csv("tissues_weight.txt", header=TRUE, check.names=FALSE, sep="\t")
ordering <- sort(tissuesWEIGHT[,"WATgon"], index.return=TRUE)$ix
sortWeight <- tissuesWEIGHT[ordering,]
plot(main="weight relationship", c(0, nrow(tissuesWEIGHT)), c(0, max(tissuesWEIGHT, na.rm=TRUE)))
lines(sortWeight[,"WATgon"], col = "red" , lwd=3 , pch=19 , type="l")
lines(sortWeight[,"Leber"], col = "blue" , lwd=3 , pch=19 , type="l")
lines(sortWeight[,"WATsc"], col = "green" , lwd=3 , pch=19 , type="l")

legend("topleft",
  legend = c("Liver", "Gon", "SCF"),
  col = c("blue", "red", "green"),
  pch = c(20,20,20),
  bty = "n",
  pt.cex = 2,
  cex = 1.2,
  text.col = "black",
  horiz = F ,
  inset = c(0.1, 0.1, 0.1))




                                                                          ## analysis on genotypes ##


genotypes <- read.csv("RawGenotypes.txt", header=TRUE, check.names=TRUE, sep="\t")
genotypes <- genotypes[,1:4]
head(genotypes)
sampleID <- unique(genotypes[,"Sample.ID"])
snpNAMES <- unique(genotypes[,"SNP.Name"])
alleles <- paste0(genotypes[,3], genotypes[,4])

for (x in 1:length(alleles)){
 if (alleles[x] == "--")
 alleles[x] = NA
}


allGenotypes <- matrix(alleles, length(snpNAMES) , length(sampleID), dimnames=list(snpNAMES, sampleID))
#allGenotypes[allGenotypes == ""] <- NA
write.table(allGenotypes, "genotypes.txt", sep = "\t", quote=FALSE)



# Dimensions we start with
dim(allGenotypes)

# QC on genotypes
# All genotypes are missing
idx <- which(apply(apply(allGenotypes, 1, is.na), 2, sum) == ncol(allGenotypes))
allGenotypes <- allGenotypes[-idx, ]
dim(allGenotypes)

# genotype is not segregating
nonSeg <- which(unlist(lapply(apply(allGenotypes,1,table), function(x){length(x) == 1})))
allGenotypes <- allGenotypes[-nonSeg, ]
dim(allGenotypes)

#quick look
mtab <- apply(allGenotypes,1,table)
mtab

# At least 2 groups with 10 observations
good <- which(unlist(lapply(apply(allGenotypes,1,table), function(x){
  length(which(x > 10)) >= 2
})))
allGenotypes <- allGenotypes[good, ]
dim(allGenotypes)

# No duplicated markers, keep the first one we see
allGenotypes <- allGenotypes[-which(duplicated(allGenotypes)),]
dim(allGenotypes)



rownames(weight) <- gsub(" ", "", rownames(weight), fixed = TRUE)


# Which individual in the phenotypes is NOT (!) in the genotype data, and throw it out of the phenotype data
indNoGeno <- which(!(rownames(weight) %in% colnames(allGenotypes)))
allGenotypes <- allGenotypes[,sort(colnames(allGenotypes))]







                                                                 ## qtl mapping ##


#bodyweight day 140
pvalsBW <- pvalsHG <- apply(allGenotypes, 1, function(geno, pheno, sex) {
  return(anova(lm(pheno ~ sex + geno))[[5]][2])
}, pheno = weight[,"140"], sex = as.numeric(id[,"Sex"]) )
plot(main = "Body weight", c(0, nrow(allGenotypes)), y = c(0, 10), t = 'n', xlab="Markers", ylab="LOD score")
points(-log10(pvalsBW))
abline(h = -log10(0.05 / nrow(allGenotypes)), col = "green")
abline(h = -log10(0.1 / nrow(allGenotypes)), col = "orange")

mAnnot <- read.csv("genotypes/SNP_Map.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE)
rownames(mAnnot) <- mAnnot[,1]
mAnnot <- mAnnot[,-c(1,4)]
mAnnot <- mAnnot[,-c(4,5,6)]
mAnnot <- mAnnot[,-3]
associatedMarkersBW <- rownames(allGenotypes)[which(-log10(pvalsBW) > 6)]
mAnnot <- mAnnot[associatedMarkersBW,]
mAnnot

#Blood glucose (Oral GTT)
pvalsOR <- apply(allGenotypes, 1, function(x){ anova(lm(oralGTT[,1] ~ x))[[5]][1] })
plot(main = "Blood glucose", c(0, nrow(allGenotypes)), y = c(0, 10), t = 'n', xlab="Markers", ylab="LOD score")
points(-log10(pvalsOR))
abline(h = -log10(0.05 / nrow(allGenotypes)), col = "green")
abline(h = -log10(0.1 / nrow(allGenotypes)), col = "orange")

mAnnot <- read.csv("genotypes/SNP_Map.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE)
rownames(mAnnot) <- mAnnot[,1]
mAnnot <- mAnnot[,-c(1,4)]
mAnnot <- mAnnot[,-c(4,5,6)]
mAnnot <- mAnnot[,-3]
associatedMarkersOR <- rownames(allGenotypes)[which(-log10(pvalsOR) > 6)]
mAnnot <- mAnnot[associatedMarkersOR,]
mAnnot


#Blood glucose (Insulin TT)
pvalsITT <- apply(allGenotypes, 1, function(x){ anova(lm(insulinTT[,4] ~ x))[[5]][1] })
plot(main = "insulin TT", c(0, nrow(allGenotypes)), y = c(0, 10), t = 'n', xlab="Markers", ylab="LOD score")
points(-log10(pvalsBW))
abline(h = -log10(0.05 / nrow(allGenotypes)), col = "green")
abline(h = -log10(0.1 / nrow(allGenotypes)), col = "orange")

mAnnot <- read.csv("genotypes/SNP_Map.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE)
rownames(mAnnot) <- mAnnot[,1]
mAnnot <- mAnnot[,-c(1,4)]
mAnnot <- mAnnot[,-c(4,5,6)]
mAnnot <- mAnnot[,-3]
associatedMarkersITT <- rownames(allGenotypes)[which(-log10(pvalsITT) > 5)]
mAnnot <- mAnnot[associatedMarkersITT,]
mAnnot

#Muscle fat
pvalsMF <- pvalsHG <- apply(allGenotypes, 1, function(geno, pheno, sex) {
  return(anova(lm(pheno ~ sex + geno))[[5]][2])
}, pheno = muscleFat[,"Fettanteil [%]"], sex = as.numeric(id[,"Sex"]) )
plot(main = "Muscle fat", c(0, nrow(allGenotypes)), y = c(0, 10), t = 'n', xlab="Markers", ylab="LOD score")
points(-log10(pvalsMF))
abline(h = -log10(0.05 / nrow(allGenotypes)), col = "green")
abline(h = -log10(0.1 / nrow(allGenotypes)), col = "orange")

mAnnot <- read.csv("genotypes/SNP_Map.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE)
rownames(mAnnot) <- mAnnot[,1]
mAnnot <- mAnnot[,-c(1,4)]
mAnnot <- mAnnot[,-c(4,5,6)]
mAnnot <- mAnnot[,-3]
associatedMarkersMF <- rownames(allGenotypes)[which(-log10(pvalsMF) > 5)]
mAnnot <- mAnnot[associatedMarkersMF,]
mAnnot


#Liver fat




#markers on chromosome 3
mAnnot <- read.csv("genotypes/SNP_Map.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE)
rownames(mAnnot) <- mAnnot[,1]
mAnnot <- mAnnot[,-c(1,4)]
mAnnot <- mAnnot[,-c(4,5,6)]
mAnnot <- mAnnot[,-3]
markerschr3 <- which(mAnnot[,1] == 3)
markerschr3 <- rownames(mAnnot[markerschr3,])
genotypesChr3 <- allGenotypes[rownames(allGenotypes) %in% markerschr3,]

#bodyweight day 140
pvalsBW <- apply(genotypesChr3, 1, function(x){ anova(lm(weight[,"140"] ~ x))[[5]][1] })
plot(c(0, nrow(genotypesChr3)), y = c(0, 10), t = 'n', xlab="markers", ylab="LOD score")
points(-log10(pvalsBW))
abline(h = -log10(0.05 / nrow(genotypesChr3)), col = "green")
abline(h = -log10(0.1 / nrow(genotypesChr3)), col = "orange")


# QTL IGTT and ORALGTT considering the area under the curve for each animals (with adjustment) #


library(DescTools)

#insulinTT
aucs <- c()
for(x in 1:nrow(insulinTT)){
  aucs <- c(aucs, AUC(x = c(0, 15, 30, 60), as.numeric(insulinTT[x, 1:4])))
}

insulinTT <- cbind(insulinTT, auc = aucs)

#calculate adjusted auc
# auc - base line auc (area from data at point 0 )

ITTaucs.adj <- cbind(apply(insulinTT[, 1:4],2, as.numeric) - mean(as.numeric(insulinTT[, 1])))

aucs.new1 <- c()
for(x in 1:nrow(ITTaucs.adj)){
  aucs.new1 <- c(aucs.new1, AUC(x = c(0, 15, 30, 60), as.numeric(ITTaucs.adj[x, 1:4])))
}

ITTaucs.adj <- cbind(ITTaucs.adj, auc = aucs.new1)


#Blood glucose (Insulin TT)
pvalsITT <- pvalsHG <- apply(allGenotypes, 1, function(geno, pheno, sex) {
  return(anova(lm(pheno ~ sex + geno))[[5]][2])
}, pheno = ITTaucs.adj[,5], sex = as.numeric(id[,"Sex"]) )
plot(main = "insulin TT", c(0, nrow(allGenotypes)), y = c(0, 10), t = 'n', xlab="Markers", ylab="LOD score")
points(-log10(pvalsITT))
abline(h = -log10(0.05 / nrow(allGenotypes)), col = "green")
abline(h = -log10(0.1 / nrow(allGenotypes)), col = "orange")

mAnnot <- read.csv("genotypes/SNP_Map.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE)
rownames(mAnnot) <- mAnnot[,1]
mAnnot <- mAnnot[,-c(1,4)]
mAnnot <- mAnnot[,-c(4,5,6)]
mAnnot <- mAnnot[,-3]
associatedMarkersITT <- rownames(allGenotypes)[which(-log10(pvalsITT) > 5)]
mAnnot <- mAnnot[associatedMarkersITT,]
mAnnot



#oralGTT
aucs <- c()
for(x in 1:nrow(oralGTT)){
  aucs <- c(aucs, AUC(x = c(0, 15, 30, 60, 120), as.numeric(oralGTT[x, 1:5])))
}

oralGTT <- cbind(oralGTT, auc = aucs)

#calculate adjusted auc
# auc - base line auc (area from data at point 0 )

oralGTTaucs.adj <- cbind(apply(oralGTT[, 1:5],2, as.numeric) - mean(as.numeric(oralGTT[, 1])))

aucs.new1 <- c()
for(x in 1:nrow(oralGTTaucs.adj)){
  aucs.new1 <- c(aucs.new1, AUC(x = c(0, 15, 30, 60, 120), as.numeric(oralGTTaucs.adj[x, 1:5])))
}

oralGTTaucs.adj <- cbind(oralGTTaucs.adj, auc = aucs.new1)

#Blood glucose (Oral GTT)
pvalsOR <- pvalsHG <- apply(allGenotypes, 1, function(geno, pheno, sex) {
  return(anova(lm(pheno ~ sex + geno))[[5]][2])
}, pheno = oralGTTaucs.adj[,6], sex = as.numeric(id[,"Sex"]) )
plot(main = "Blood glucose", c(0, nrow(allGenotypes)), y = c(0, 10), t = 'n', xlab="Markers", ylab="LOD score")
points(-log10(pvalsOR))
abline(h = -log10(0.05 / nrow(allGenotypes)), col = "green")
abline(h = -log10(0.1 / nrow(allGenotypes)), col = "orange")

mAnnot <- read.csv("genotypes/SNP_Map.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE)
rownames(mAnnot) <- mAnnot[,1]
mAnnot <- mAnnot[,-c(1,4)]
mAnnot <- mAnnot[,-c(4,5,6)]
mAnnot <- mAnnot[,-3]
associatedMarkersOR <- rownames(allGenotypes)[which(-log10(pvalsOR) > 4)]
mAnnot <- mAnnot[associatedMarkersOR,]
mAnnot












##manhattanPlot (not working yet)
mAnnot <- read.csv("genotypes/SNP_Map.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE)
rownames(mAnnot) <- mAnnot[,1]
mAnnot <- mAnnot[,-c(1,4)]
mAnnot <- mAnnot[,-c(4,5,6)]
mAnnot <- mAnnot[,-3]
pvalsOR <- data.frame(pvalsOR)
markers <- rownames(pvalsOR)
markersAnnot <- mAnnot[markers,]
manhattanDataset <- cbind(markersAnnot, pvalsOR)
colnames(manhattanDataset) <- c("Chr", "Pos", "marPvalue")
manhattanDataset <- manhattanDataset[order("Chr"),]

#manhattanDataset, contains columns named: Chr, Pos, marPvalue
plotLODcurves <- function(manhattanDataset, gap = 50000000, chromosomes = c(1:19, "X")){
  manhattanSorted <- NULL
  chr.starts <- c(0)
  prev <- 0

  for(chr in chromosomes){     #defining the start of each chromosome
    onchr <- manhattanDataset[manhattanDataset$Chr == chr,]
    manhattanSorted <- rbind(manhattanSorted, onchr)
    chr.starts <- c(chr.starts, gap + prev + max(as.numeric(as.character(onchr$Pos))))
    prev <- chr.starts[length(chr.starts)]
  }

  plot(c(0, chr.starts[length(chr.starts)]), y = c(min(onchr$marPvalue), max(onchr$marPvalue) * 5), xaxs="i", yaxs="i", t = 'n', xaxt="n", xlab="Chromsome", ylab="LOD score")
  i <- 1
  abline(h = seq(5, 50, 5), lty=2, lwd=1, col=rgb(0.5, 0.5, 0.5, 0.5))
  mCol <- c("black", "orange")
  for(chr in chromosomes){
    onchr <- manhattanSorted[manhattanSorted$Chr == chr,]
    points(x = onchr$Pos + chr.starts[i], y = onchr$marPvalue, t = "p", pch = 1, col=mCol[as.numeric((i %%2) == 0)+1])
    i <- i + 1
  }
  axis(1, at = (chr.starts  + 0.5 * diff(chr.starts))[1:length(chromosomes)], chromosomes)

  for (x in 1:length(chr.starts))
  abline(v=chr.starts[x])
  for (n in 2:length(chr.starts))
  abline(v=chr.starts[n] - gap)
}



op <- par(mfrow=c(2,1))
plotLODcurves(manhattanDataset)
manhattanDataset$marPvalue <- rnorm(length(manhattanDataset$marPvalue), 5, 100)
plotLODcurves(manhattanDataset)
