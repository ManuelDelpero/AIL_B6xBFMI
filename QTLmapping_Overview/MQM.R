# AIL_S1xS2 Analysis on Phenotypes
#
# copyright (c) - Manuel Delpero
# first written september, 2019
# modelling using AIL_B6xS1

setwd("C:/Users/Manuel/Desktop/AIL_B6xBFMI/RAWDATA")

genotypes <- read.csv("genomatrix.clean_numeric.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
phenotypes <- read.csv("allPhenotypes.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
markerannot <- read.csv("SNP_Map.txt", header=TRUE, sep="\t", row.names=2, check.names=FALSE)
markerannot <- markerannot[,-1]
markerannot <- markerannot[rownames(genotypes),]
markerannot <- markerannot[sort(markerannot[,"Position"], index.return=TRUE)$ix,]

chromosomes <- c(1:19, "X", "Y")

annotation <- c()
for(chr in chromosomes){
  annotation <- rbind(annotation, markerannot[markerannot[,"Chromosome"] == chr,])
}

phenonames <- colnames(phenotypes)[-c(1:5)]

# Make sure that the ordering between phenotypes and genotypes matches !!!!!
# Also sort the markers by their natural chromosome ordering
genotypes <- genotypes[rownames(annotation), rownames(phenotypes)]
write.table(genotypes, "OrderedGenotypes.txt", sep = "\t", quote=FALSE)

# Covariates we could/need to include in the model, we test them on their pvalue
sex <- phenotypes[, "Sex"]
wg <- phenotypes[, "WG"]
mother <- phenotypes[, "Mutter"]

#MQM for day 112 bodyweight
d112 <- apply(genotypes, 1, function(geno, pheno, sex, markerChr3, mother){
  return(anova(lm(pheno ~ sex + markerChr3 + mother + geno))["Pr(>F)"]["geno",])
  }, pheno = phenotypes[,"d112"], sex = sex, markerChr3 = as.factor(genotypes["gUNC5036315",]) , mother = mother)

lodscores <- -log10(d112)
plot(lodscores, col = as.numeric(as.factor(annotation[,"Chromosome"])))
abline(h = -log10(0.05 / length(lodscores)), col="orange")
abline(h = -log10(0.01 / length(lodscores)), col="green")

# QQplot
pvals.exp <- (rank(d112, ties.method="first")+0.5) / (length(d112) + 1)
plot(-log10(pvals.exp), -log10(d112))
abline(a=0, b=1)