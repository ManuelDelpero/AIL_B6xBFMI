setwd("C:/Github/AIL_B6xBFMI/QTLmapping_GlucoseTests")

library(DescTools)
genotypes <- read.csv("genomatrix.clean.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
phenotypes <- read.csv("phenotypes.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
insulinTT <- phenotypes[, ]
oralGTT <- phenotypes[, ]

# QTL mapping using data from insulinTT
# calculate the area under the curve (auc) for each animal
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
pvalsITT <- pvalsHG <- apply(genotypes, 1, function(geno, pheno, sex) {
  return(anova(lm(pheno ~ sex + geno))[[5]][2])
}, pheno = ITTaucs.adj[,5], sex = as.numeric(phenotypes[,"Sex"]) )
plot(main = "insulin TT", c(0, nrow(allGenotypes)), y = c(0, 10), t = 'n', xlab="Markers", ylab="LOD score")
points(-log10(pvalsITT))
abline(h = -log10(0.05 / nrow(allGenotypes)), col = "green")
abline(h = -log10(0.1 / nrow(allGenotypes)), col = "orange")

## QTL mapping using data from oralGTT
# calculate the area under the curve (auc) for each animal
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
}, pheno = oralGTTaucs.adj[,6], sex = as.numeric(phenotypes[,"Sex"]) )
plot(main = "Blood glucose", c(0, nrow(allGenotypes)), y = c(0, 10), t = 'n', xlab="Markers", ylab="LOD score")
points(-log10(pvalsOR))
abline(h = -log10(0.05 / nrow(allGenotypes)), col = "green")
abline(h = -log10(0.1 / nrow(allGenotypes)), col = "orange")
