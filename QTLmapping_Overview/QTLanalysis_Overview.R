setwd("D:/Edrive/Mouse/AIL_Manuel")

genotypes <- read.csv("genomatrix.clean.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
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

phenonames <- colnames(phenotypes)[-c(1:3)]

# Make sure that the ordering between phenotypes and genotypes matches !!!!!
# Also sort the markers by their natural chromosome ordering
genotypes <- genotypes[rownames(annotation), rownames(phenotypes)]

pmatrix <- matrix(NA, nrow(genotypes), length(phenonames), dimnames= list(rownames(genotypes), phenonames))
shapmatrix <- matrix(NA, nrow(genotypes), length(phenonames), dimnames= list(rownames(genotypes), phenonames))
for (pname in phenonames) {
  pvalues <- apply(genotypes, 1, function(geno, pheno, sex) {
    mmodel <- lm(pheno ~ geno)
    shaptest <- shapiro.test(mmodel$residuals)
    return(c(anova(mmodel)[[5]][1], shaptest$p.value))
  }, pheno = phenotypes[,pname], sex = phenotypes[,"Sex"])
  pmatrix[colnames(pvalues), pname] <- pvalues[1,]
  shapmatrix[colnames(pvalues), pname] <- pvalues[2,]
  cat("Done: ", pname, "\n")
}
lodmatrix <- -log10(pmatrix)

write.table(lodmatrix, "lodmatrix.txt", sep = "\t", quote=FALSE)

for (pname in phenonames) {
  plot(lodmatrix[,pname], main = pname, col = as.numeric(as.factor(annotation[,"Chromosome"])))
  abline(h = -log10(0.05 / nrow(lodmatrix)), col="orange")
  abline(h = -log10(0.01 / nrow(lodmatrix)), col="green")
}

signmatrix <- lodmatrix[which(apply(lodmatrix, 1, function(x){ any(x > -log10(0.01 / nrow(lodmatrix))) })),]
signannotmatrix <- cbind(annotation[rownames(signmatrix), ], signmatrix)

write.table(signannotmatrix, "signannotmatrix.txt", sep = "\t", quote=FALSE)
