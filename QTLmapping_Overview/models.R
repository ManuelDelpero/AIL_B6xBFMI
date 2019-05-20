setwd("C:/Users/Manuel/Desktop/AIL_B6xBFMI/RAWDATA")

genotypes <- read.csv("OrderedGenotypes.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
phenotypes <- read.csv("allPhenotypes.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
annotation <- read.csv("annotation.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE)

# correlation between the phenotypes
phenoCor <- cor(data.matrix(phenotypes), use="pairwise.complete.obs")
write.table(phenoCor, "CorrelationPhenotypes.txt", sep = "\t", quote=FALSE)

# best model for triglycerides 
anova(lm(phenotypes[,"Triglycerides/Proteins"] ~ phenotypes[,"WATsc"]))
anova(lm(phenotypes[,"Triglycerides/Proteins"] ~ phenotypes[,"d140"]))

Triglicerydes <- apply(genotypes,1,function(geno, pheno){
 return(anova(lm(pheno ~ geno))[[5]][1])
}, pheno = phenotypes[,"Triglycerides/Proteins"])

#Triglicerydes <- apply(genotypes,1,function(geno, pheno, phenocor1, phenocor2){
 #return(anova(lm(pheno ~ phenocor1 + phenocor2 + geno))[[5]][1])
#}, pheno = phenotypes[,"Triglycerides/Proteins"], phenocor1 = phenotypes[,"WATsc"], phenocor2 = phenotypes[,"d140"])

lodscores <- -log10(Triglicerydes)
plot(lodscores, col = as.numeric(as.factor(annotation[,"Chromosome"])))
abline(h = -log10(0.05 / length(lodscores)), col="orange")

# QQplot
pvals.exp <- (rank(Triglicerydes, ties.method="first")+0.5) / (length(Triglicerydes) + 1)
plot(-log10(pvals.exp), -log10(Triglicerydes))
abline(a=0, b=1)

# Lambda correction
lambda.B <- round(median(qchisq(1.0 - Triglicerydes, 1), na.rm=TRUE) /  qchisq(0.5, 1),3)
lambda.B 
# no correction needed
signscores <- lodscores[which(lodscores > -log10(0.05 / length(p.adj)))]
signscores <- lodscores[which(lodscores > -log10(0.05 / length(p.adj)))]
signannot <- annotation[names(signscores),]


# best model for liver weight
anova(lm(phenotypes[,"Leber"] ~ phenotypes[,"Length"]))
anova(lm(phenotypes[,"Leber"] ~ phenotypes[,"d140"]))

LiverWeight <- apply(genotypes,1,function(geno, pheno, sex, mother, Length){
 return(anova(lm(pheno ~ sex + mother + Length + geno))[[5]][4])
}, pheno = phenotypes[,"Leber"], mother = phenotypes[,"Mutter"], sex = phenotypes[,"Sex"], Length = phenotypes[,"Length"])

lodscores <- -log10(LiverWeight)
plot(lodscores, col = as.numeric(as.factor(annotation[,"Chromosome"])))
abline(h = -log10(0.05 / length(lodscores)), col="orange")

# QQplot
pvals.exp <- (rank(LiverWeight, ties.method="first")+0.5) / (length(LiverWeight) + 1)
plot(-log10(pvals.exp), -log10(LiverWeight))
abline(a=0, b=1)

# Lambda correction
lambda.B <- round(median(qchisq(1.0 - LiverWeight, 1), na.rm=TRUE) /  qchisq(0.5, 1),3)
# Adjusted pvalues
chiSq <- qchisq(LiverWeight, 1, lower.tail = FALSE)
LiverWeight <- pchisq(chiSq / lambda.B, 1, lower.tail = FALSE)
plot(-log10(LiverWeight ), col = as.numeric(as.factor(annotation[,"Chromosome"])))
abline(h = -log10(0.05 / length(LiverWeight )), col="orange")
lambda.A <- round(median(qchisq(1.0 - LiverWeight, 1), na.rm=TRUE) /  qchisq(0.5, 1),3)
# QQplot adj 
pvals.exp <- (rank(LiverWeight , ties.method="first")+0.5) / (length(LiverWeight ) + 1)
plot(-log10(pvals.exp), -log10(LiverWeight))
abline(a=0, b=1)
lodscores <- -log10(LiverWeight)
signscores <- lodscores[which(lodscores > -log10(0.05 / length(LiverWeight )))]
signannot <- annotation[names(signscores),]


