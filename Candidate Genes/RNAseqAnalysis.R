# AIL_S1xS2 analysis of RNA-seq results, check if the transcript ENSMUST00000125471 is expressed in BFMI and B6 in liver
#
# copyright (c) - Manuel Delpero
# first written july, 2020
# 

setwd("C:/Users/Manuel/Desktop/AIL_B6xBFMI/RAWDATA")

TransExpression <- read.table("transcript_expression.txt", sep = "\t", header = TRUE)
annotation <- read.table("SampleDescription.txt", sep = "\t", check.names = FALSE, header = TRUE)

annotation <- annotation[c(14,16),]
colnames(TransExpression) <- gsub(".aln.sort.rmdup.rg.realigned.recal.bam", "",colnames(TransExpression))
colnames(TransExpression) <- gsub("X", "",colnames(TransExpression))

ExpressionLiverBFMI <- TransExpression[,c(7,8)]
ExpressionLiverB6 <- TransExpression[,c(9,10)]
ExpressionGonBFMIB6 <- TransExpression[,c(27,29)]

# Check the expression of the transcript ENSMUST00000125471 of FTO (Stop gained mutation in our BFMI lines)

ExpressionLiverBFMI[which(rownames(ExpressionLiverBFMI) == "ENSMUST00000125471"),] # it is expressed in BFMI!
ExpressionLiverB6[which(rownames(ExpressionLiverB6) == "ENSMUST00000125471"),] # it is expressed in B6 but seems downregulated in males but not females
ExpressionGonBFMIB6[which(rownames(ExpressionGonBFMIB6) == "ENSMUST00000125471"),]


