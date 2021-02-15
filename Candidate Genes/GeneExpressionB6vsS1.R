setwd("C:/Users/Manuel/Desktop/AIL_B6xBFMI/RAWDATA/Microarrays_S1")
arraymapping <- read.table("mapping.txt", sep = '\t', header=TRUE, colClasses = "character", row.names = 1)

# Process S1 arrays
library(affy)
library("gplots")

dat <- ReadAffy(cdfname ='ClariomSMouse_Mm_ENST') # Use the clariomsmouse CDF from http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/23.0.0/enst.asp
eset <- mas5(dat)

library("affyPLM")

# QC 
expressions <- log2(assayData(eset)$exprs)
expressions <- normalize.quantiles(expressions)
rownames(expressions) <- rownames(assayData(eset)$exprs)
colnames(expressions) <- colnames(assayData(eset)$exprs)

ids <- unlist(lapply(strsplit(colnames(expressions), "_"),"[",2))

liver <- which(grepl("L", ids))
gonadalfat <- which(grepl("G", ids))
skeletalmuscle <- which(grepl("S", ids))
pankreas <- which(grepl("P", ids))

colnames(expressions) <- arraymapping[ids, 1]

corM <- cor(expressions)
image(corM)

# Get only the exression for liver
Liver <- grep("L", colnames(expressions))
expressionsLS1 <- expressions[,Liver]


# Process B6 arrays
setwd("C:/Users/Manuel/Desktop/AIL_B6xBFMI/RAWDATA/Microarrays_B6")
arraymapping <- read.table("mapping.txt", sep = '\t', header=TRUE, colClasses = "character")


library(affy)
library("gplots")

dat <- ReadAffy(cdfname ='ClariomSMouseHT_Mm_ENST')
eset <- mas5(dat)

library("affyPLM")

# QC 
expressions <- log2(assayData(eset)$exprs)
expressions <- normalize.quantiles(expressions)
rownames(expressions) <- rownames(assayData(eset)$exprs)
colnames(expressions) <- colnames(assayData(eset)$exprs)

corM <- cor(expressions)
expressionsLB6 <- expressions

# Combine the two datasets 
ExpressionsS1B6L <- cbind(expressionsLS1[,c(1:7)], expressionsLB6[,c(1:8)])

# Some quality check
corM <- cor(ExpressionsS1B6L)
heatmap(corM) 				# one outlier in S1 liver samples

colnames(ExpressionsS1B6L) <- c("S1-7247L", "S1-7248L", "S1-7249L", "S1-7298L", "S1-7299L", "S1-7304L", "S1-7305L", "1-B6L","2-B6L", "3-B6L", "4-B6L", "5-B6L", "6-B6L", "7-B6L", "8-B6L")
corM <- cor(ExpressionsS1B6L)
heatmap(corM) 		

# Diff. gene expression analysis
getSignificant <- function(expressions, Tissue = "L", adjust = "BH", p.val = 0.05){
  S1P <- which(grepl("S1", colnames(expressions)) & grepl(Tissue, colnames(expressions)))
  B6P <- which(grepl("B6", colnames(expressions)) & grepl(Tissue, colnames(expressions)))

  res <- t(apply(expressions[, c(S1P, B6P)],1, function(x) {
    s1v <- x[1:length(S1P)]
    B6v <- x[(length(S1P)+1):(length(S1P) + length(B6P))]
    pval <- tryCatch({t.test(s1v, B6v)$p.value}, error = function(e) { return(NA) })
    return(c(mean(s1v), mean(B6v), sd(s1v), sd(B6v), mean(s1v) / mean(B6v), log2(mean(s1v) / mean(B6v)), pval))
  }))
  colnames(res) <- c("mean(s1)", "mean(B6)", "sd(s1)", "sd(s2)", "FC", "logFC", "p.value")
  significant <- res[which(p.adjust(res[,"p.value"], adjust) < p.val),]
  rownames(significant) <- gsub("_at", "", rownames(significant))
  return(significant)
}

annotate <- function(significant){
  library(biomaRt)
  bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

  res.biomart <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id", "mgi_id", "mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position", "strand"), 
                       filters = c("ensembl_transcript_id_version"), 
                       values = rownames(significant), mart = bio.mart)
  rownames(res.biomart) <- res.biomart[,1]
  annotated <- cbind(res.biomart[rownames(significant), ], significant)
  return(annotated)
}

liverExpression <- annotate(getSignificant(expressionS1B6L, "L", p.val = 0.05)) 		# Bbs7 is downregulated in S1 
