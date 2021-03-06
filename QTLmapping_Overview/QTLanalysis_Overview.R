setwd("C:/Users/Manuel/Desktop/AIL_B6xBFMI/RAWDATA")
# for QTL mapping with forward and backward model selection
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin Manuel Delpero & Danny Arends
# 
# first written november, 2019

#genotypes <- read.csv("genomatrix.clean.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
genotypes <- read.csv("genomatrix.clean_numeric.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
phenotypes <- read.csv("PhenotypesComplete.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
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

# Covariates we could/need to include in the model, we test them on their pvalue, if significant they are added into the model, forward covariate selection
sex <- phenotypes[, "Sex"]
wg <- phenotypes[, "WG"]
mother <- phenotypes[, "Mutter"]

pmatrix <- matrix(NA, nrow(genotypes), length(phenonames), dimnames= list(rownames(genotypes), phenonames))
for (pname in phenonames){
  pheno <- phenotypes[, pname]
  p.sex <- anova(lm(pheno ~ sex))["Pr(>F)"]["sex",]
  p.mother <- anova(lm(pheno ~ mother))["Pr(>F)"]["mother",]
  p.wg <- anova(lm(pheno ~ wg))["Pr(>F)"]["wg",]
  myfactors <- c()
  if(p.sex < 0.05) myfactors <- c(myfactors, "sex")
  if(p.mother < 0.05) myfactors <- c(myfactors, "mother")
  if(p.wg < 0.05) myfactors <- c(myfactors, "wg")
  myformula <- paste0("pheno ~ ", paste(c(myfactors, "geno"), collapse = " + "))
  cat(pname, " ", myformula, "\n")
  
  pvalues <- apply(genotypes, 1, function(geno, pheno, sex, wg, mother, myformula) {
    mmodel <- lm(formula(myformula))
    return(anova(mmodel)["Pr(>F)"]["geno",])
  }, pheno = phenotypes[,pname], sex = sex, wg = wg, mother = mother, myformula = myformula)
  pmatrix[names(pvalues), pname] <- pvalues
}
lodmatrix <- -log10(pmatrix)

write.table(lodmatrix, "lodmatrix.txt", sep = "\t", quote=FALSE)

# Manhattan plots
for (pname in  phenonames){   #phenonames) {
  plot(lodmatrix[,pname], main = pname, col = as.numeric(as.factor(annotation[,"Chromosome"])))
  abline(h = -log10(0.05 / nrow(lodmatrix)), col="orange")
  abline(h = -log10(0.01 / nrow(lodmatrix)), col="green")
}

# QQplots
for (pname in "Plasma_Triglycerides"){   #phenonames) {
  pvals.exp <- (rank(pmatrix[,pname], ties.method="first")+0.5) / (length(pmatrix[,pname]) + 1)
  plot(-log10(pvals.exp), -log10(pmatrix[,pname]))
  abline(a=0, b=1)
}

# Backward covariate selection to include only the covariates that actually improve the model according to the QQplots,
# and that still have a significant influence when tested in the full model

#phenonames <- c("WATgon", "Leber", "WATsc", "Plasma_Trigylcerides", "Plasma_Cholesterol")
anova(lm(phenotypes[,"WATgon"] ~ phenotypes[,"Sex"] + phenotypes[,"WG"] )) # exclude WG
mmodelWATgon <- lm(phenotypes[,"WATgon"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
mmodelWATsc <- lm(phenotypes[,"WATsc"] ~ phenotypes[,"Mutter"])
mmodelLeber <- lm(phenotypes[,"Leber"] ~ phenotypes[,"Sex"] + phenotypes[,"Mutter"])
mmodelPlasmaTrig <- lm(phenotypes[,"Plasma_Triglycerides"] ~ phenotypes[,"Sex"])
anova(lm(phenotypes[,"Plasma_Cholesterol"] ~ phenotypes[,"Sex"] +  phenotypes[,"Mutter"])) # try to exclude the mother and perform lambda correction

pmatrix <- matrix(NA, nrow(genotypes), length(phenonames), dimnames= list(rownames(genotypes), phenonames))
for (pname in phenonames){
  pheno <- phenotypes[, pname]
  p.sex <- anova(lm(pheno ~ sex))["Pr(>F)"]["sex",]
  p.mother <- anova(lm(pheno ~ mother))["Pr(>F)"]["mother",]
  p.wg <- anova(lm(pheno ~ wg))["Pr(>F)"]["wg",]
  myfactors <- c()
  if(p.sex < 0.05) myfactors <- c(myfactors, "sex")
  if((p.mother < 0.05) && (pname != "ITTauc") && (pname != "WATgon")) myfactors <- c(myfactors, "mother")
  if((p.wg < 0.05) && (pname != "WATgon")) myfactors <- c(myfactors, "wg")
  myformula <- paste0("pheno ~ ", paste(c(myfactors, "geno"), collapse = " + "))
  cat(pname, " ", myformula, "\n")
  
  pvalues <- apply(genotypes, 1, function(geno, pheno, sex, wg, mother, myformula) {
    mmodel <- lm(formula(myformula))
    return(anova(mmodel)["Pr(>F)"]["geno",])
  }, pheno = phenotypes[,pname], sex = sex, wg = wg, mother = mother, myformula = myformula)
  pmatrix[names(pvalues), pname] <- pvalues
}
lodmatrix <- -log10(pmatrix)

# Have a look at the results with the new models
# Manhattan plots
for (pname in  "ITTauc"){   #phenonames) {
  plot(lodmatrix[,pname], main = pname, col = as.numeric(as.factor(annotation[,"Chromosome"])))
  abline(h = -log10(0.05 / nrow(lodmatrix)), col="orange")
  abline(h = -log10(0.01 / nrow(lodmatrix)), col="green")
}

# QQplots
for (pname in "ITTauc"){   #phenonames) {
  pvals.exp <- (rank(pmatrix[,pname], ties.method="first")+0.5) / (length(pmatrix[,pname]) + 1)
  plot(-log10(pvals.exp), -log10(pmatrix[,pname]))
  abline(a=0, b=1)
}

signmatrix <- lodmatrix[which(apply(lodmatrix, 1, function(x){ any(x > -log10(0.05 / nrow(lodmatrix))) })),]
signannotmatrix <- cbind(annotation[rownames(signmatrix), ], signmatrix)

write.table(signannotmatrix, "signannotmatrix.txt", sep = "\t", quote=FALSE)

# lambda correction for each phenotype
pvalues.adj <- matrix(NA, nrow(genotypes), length(phenonames), dimnames= list(rownames(genotypes), phenonames))
for (pheno in colnames(pvalues.adj)){
  lambda.B <- round(median(qchisq(1.0 - pmatrix[,pheno], 1), na.rm=TRUE) /  qchisq(0.5, 1),3)
    if (lambda.B > 1.15){
      chiSq <- qchisq(pmatrix[,pheno], 1, lower.tail = FALSE)
      p.adj <- pchisq(chiSq / lambda.B, 1, lower.tail = FALSE)
    }else{ p.adj <- pmatrix[,pheno] }
      pvalues.adj[names(p.adj), pheno] <- p.adj
}
lodmatrix.adj <- -log10(pvalues.adj)

# additional models	for liver weight and Triglycerides
source("C:/Github/AIL_B6xBFMI/QTLmapping_Overview/models.R")
# creating the pvalues matrix with adjusted pvalues
pmatrix.adj <- cbind(pvalues.adj, LiverWeight)
lodmatrix.adj <- -log10(pmatrix.adj)
write.table(lodmatrix.adj, "lodmatrix.adj.txt", sep = "\t", quote=FALSE)
signmatrix.adj <- lodmatrix.adj[which(apply(lodmatrix.adj, 1, function(x){ any(x > -log10(0.05 / nrow(lodmatrix.adj))) })),]
signannotmatrix.adj <- cbind(annotation[rownames(signmatrix.adj), ], signmatrix.adj)
write.table(signannotmatrix.adj, "signannotmatrix.adj.txt", sep = "\t", quote=FALSE)

getVarianceExplained <- function(genotypes, phenotypes, pheno.col = "d77", marker = "gUNC5036315"){
  genotype <- as.factor(t(genotypes[marker,]))
  littersize <- as.factor(phenotypes[, "WG"])
  subfamily <- as.factor(phenotypes[, "Mutter"])
  sex <- as.factor(phenotypes[, "Sex"])
  phenotype <- phenotypes[, pheno.col]
  model <- lm(phenotype ~ subfamily + littersize + sex + genotype)
  tryCatch(res  <- anova(model), error = function(e){ res <<- NA })
  cat(names(model$coefficients),"\n")
  cat(model$coefficients,"\n")
  varExplained  <- res[, "Sum Sq"] / sum((phenotype - mean(phenotype, na.rm=TRUE))^2, na.rm=TRUE)
  names(varExplained) <- c("subfamily","littersize","sex", "marker")
  return(round(varExplained * 100, digits=1))
}

getVarianceExplained(genotypes, phenotypes, pheno.col = "Triglycerides/Proteins", marker = "S1H083826428") # 20.5 % of the variance explained by the QTL on chr 8!!