setwd("D:/Edrive/Mouse/AIL_Manuel")

# read files

## analysis on genotypes ##

## qtl mapping ##

#bodyweight day 140
pvalsBW <- pvalsHG <- apply(allGenotypes, 1, function(geno, pheno, sex) {
  return(anova(lm(pheno ~ sex + geno))[[5]][2])
}, pheno = weight[,"140"], sex = as.numeric(id[,"Sex"]) )
plot(main = "Body weight", c(0, nrow(allGenotypes)), y = c(0, 10), t = 'n', xlab="Markers", ylab="LOD score")
points(-log10(pvalsBW))
abline(h = -log10(0.05 / nrow(allGenotypes)), col = "green")
abline(h = -log10(0.1 / nrow(allGenotypes)), col = "orange")

mAnnot <- read.csv("SNP_Map.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE)
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

mAnnot <- read.csv("SNP_Map.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE)
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

mAnnot <- read.csv("SNP_Map.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE)
rownames(mAnnot) <- mAnnot[,1]
mAnnot <- mAnnot[,-c(1,4)]
mAnnot <- mAnnot[,-c(4,5,6)]
mAnnot <- mAnnot[,-3]
associatedMarkersITT <- rownames(allGenotypes)[which(-log10(pvalsITT) > 5)]
mAnnot <- mAnnot[associatedMarkersITT,]
mAnnot

#markers on chromosome 3
mAnnot <- read.csv("SNP_Map.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE)
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
