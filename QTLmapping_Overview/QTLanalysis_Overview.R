setwd("C:/Users/Manuel/Desktop/AIL_B6xBFMI/RAWDATA")

#genotypes <- read.csv("genomatrix.clean.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
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

pmatrix <- matrix(NA, nrow(genotypes), length(phenonames), dimnames= list(rownames(genotypes), phenonames))
shapmatrix <- matrix(NA, nrow(genotypes), length(phenonames), dimnames= list(rownames(genotypes), phenonames))
for (pname in phenonames) {
  pvalues <- apply(genotypes, 1, function(geno, pheno, sex, wg) {
    mmodel <- lm(pheno ~ sex + geno + wg)
    shaptest <- shapiro.test(mmodel$residuals)
    return(c(anova(mmodel)["Pr(>F)"]["geno",], shaptest$p.value))
  }, pheno = phenotypes[,pname], sex = phenotypes[,"Sex"], wg = phenotypes[,"WG"])
  pmatrix[colnames(pvalues), pname] <- pvalues[1,]
  shapmatrix[colnames(pvalues), pname] <- pvalues[2,]
  cat("Done: ", pname, "\n")
}
lodmatrix <- -log10(pmatrix)

write.table(lodmatrix, "lodmatrix.txt", sep = "\t", quote=FALSE)

#Manhattan plots
for (pname in phenonames) {
  plot(lodmatrix[,pname], main = pname, col = as.numeric(as.factor(annotation[,"Chromosome"])))
  abline(h = -log10(0.05 / nrow(lodmatrix)), col="orange")
  abline(h = -log10(0.01 / nrow(lodmatrix)), col="green")
}

#QQplots
for (pname in phenonames)
 pvals.exp <- (rank(lodmatrix[,pname], ties.method="first")+0.5) / (length(lodmatrix[,pname]) + 1)
 plot(-log10(pvals.exp), -log10(lodmatrix[,pname]))
 abline(a=0, b=1)

signmatrix <- lodmatrix[which(apply(lodmatrix, 1, function(x){ any(x > -log10(0.01 / nrow(lodmatrix))) })),]
signannotmatrix <- cbind(annotation[rownames(signmatrix), ], signmatrix)

write.table(signannotmatrix, "signannotmatrix.txt", sep = "\t", quote=FALSE)

# Analysis on the Top markers
# read the normal genotypes (non numeric)
genotypes <- read.csv("genomatrix.clean.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")

# ITT TOP marker
topmarker <- t(genotypes["gUNCHS024558",])
genophenoITT <- cbind(topmarker, phenotypes[,"ITT_0"])
colnames(genophenoITT) <- c("Genotype", "ITT")
boxplot(as.numeric(as.character(genophenoITT[, "ITT"]))  ~ genophenoITT[,"Genotype"], main = "ITT" , xlab = "Genotype", ylab = "AUC (adj)", col = (c("gold" , "darkgreen" , "lightblue")))

# Body weight TOP marker
topmarker <- t(genotypes["gUNC5046545",])
genophenoWeight <- cbind(topmarker, phenotypes[,"d98"])
colnames(genophenoWeight) <- c("Genotype", "d98")
boxplot(as.numeric(as.character(genophenoWeight[, "d98"]))  ~ genophenoWeight[,"Genotype"],  main = "Body weight" , xlab = "Genotype", ylab = "Weight (grams)", col = (c("gold" , "darkgreen" , "lightblue")))

# Quadriceps weight TOP marker
topmarker <- t(genotypes["SAH033394051",])
genophenoQuadri <- cbind(topmarker, phenotypes[,"Quadri"])
colnames(genophenoQuadri) <- c("Genotype", "Quadri")
boxplot(as.numeric(as.character(genophenoQuadri[, "Quadri"]))  ~ genophenoQuadri[,"Genotype"],  main = "Quadriceps weight" , xlab = "Genotype", ylab = "Weight (grams)", col = (c("gold" , "darkgreen")))
