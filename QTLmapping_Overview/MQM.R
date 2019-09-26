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

# Effect plot with normal genotypes
genotypes <- read.csv("genomatrix.clean.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
topmarker <- t(genotypes["gUNC10595065",])
genophenoWeight <- cbind(topmarker, phenotypes[,"d112"])
colnames(genophenoWeight) <- c("Genotype", "d112")
boxplot(as.numeric(as.character(genophenoWeight[, "d112"]))  ~ genophenoWeight[,"Genotype"],  main = "Body weight" , xlab = "Genotype", ylab = "Weight (grams)", col = (c("gold" , "darkgreen" , "lightblue")))

# LOD curve
lodmatrix <- data.frame(lodscores)
markerannot <- markerannot[rownames(genotypes),]
markerannot <- markerannot[sort(markerannot[,"Position"], index.return=TRUE)$ix,]
lodannotmatrix <- cbind(annotation[rownames(lodmatrix), ], lodmatrix)
dataset <- lodannotmatrix[, c("Chromosome", "Position","lodscores")]
chr6 <- dataset[which(dataset[,"Chromosome"] == 6),]

plot(main = "Lod score curve Chr 6", c(min(as.numeric(chr6[, "Position"])), max(as.numeric(chr6[, "Position"]))), c(0,9), ylab = "-log10(pvalue)", xlab = "Position", las = 2, t = "n", xaxt = "n")
  points(x = as.numeric(chr6[,"Position"]), y = chr6[,"lodscores"] , type = "l", col="dodgerblue", lwd = 1)
  abline(h=4.5, col="green")
  abline(h=4, col="orange")
  axis(1, at = c(0,25000000, 50000000, 75000000, 100000000, 125000000, 150000000), c("0 mb", "25 mb", "50 mb", "75 mb", "100 mb", "125 mb", "150 mb"))