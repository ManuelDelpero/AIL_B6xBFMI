setwd("C:/Users/Manuel/Desktop/AIL_B6xBFMI/RAWDATA/Microarrays")
arraymapping <- read.table("mapping.txt", sep = '\t', header=TRUE, colClasses = "character", row.names=1)

library("gplots")
library(affy)

# dat <- ReadAffy(cdfname ='ClariomSMouse_Mm_ENST') # Use the clariomsmouse CDF from http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/23.0.0/enst.asp
dat <- ReadAffy(cdfname ='clariomsmousemmenstcdf')
eset <- mas5(dat)

library("affyPLM")

# QC 
expressions <- log2(assayData(eset)$exprs)
expressions <- normalize.quantiles(expressions)
rownames(expressions) <- rownames(assayData(eset)$exprs)
colnames(expressions) <- colnames(assayData(eset)$exprs)

ids <- unlist(lapply(strsplit(colnames(expressions), "_"),"[",2))
ids <- ids[-(58:65)]
liver <- which(grepl("L", ids))
gonadalfat <- which(grepl("G", ids))
skeletalmuscle <- which(grepl("S", ids))
pankreas <- which(grepl("P", ids))

colnames(expressions) <- c(arraymapping[ids, 1], "B6-MET1L","B6-MET2L","B6-MET3L","B6-MET4L","B6-MET5L","B6-MET6L","B6-1L", "B6-2L")
corM <- cor(expressions)

heatmap(corM)

# Select the liver samples in B6 and S1 and perform the analysis in this samples (Use the results for prioritization of the genes in the regions)
expressions <- expressions[,c(1:7,58:65)]
expressions <- expressions[,-c(2, 11:13)]
corM <- cor(expressions)

heatmap(corM)

# Diff. gene expression analysis
getSignificant <- function(expressions, Tissue = "L", adjust = "BH", p.val = 0.05){
  S1P <- which(grepl("S1", colnames(expressions)) & grepl(Tissue, colnames(expressions)))
  B6P <- which(grepl("B6", colnames(expressions)) & grepl(Tissue, colnames(expressions)))

  res <- t(apply(expressions[, c(S1P, B6P)],1, function(x) {
    s1v <- x[1:length(S1P)]
    s2v <- x[(length(S1P)+1):(length(S1P) + length(B6P))]
    pval <- tryCatch({t.test(s1v, s2v)$p.value}, error = function(e) { return(NA) })
    return(c(mean(s1v), mean(s2v), sd(s1v), sd(s2v), mean(s1v) / mean(s2v), log2(mean(s1v) / mean(s2v)), pval))
  }))
  colnames(res) <- c("mean(s1)", "mean(s2)", "sd(s1)", "sd(s2)", "FC", "logFC", "p.value")
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

liver <- annotate(getSignificant(expressions, "L", p.val = 0.05))

write.table(liver, "liver_expressions_ann.txt", sep = "\t", quote = FALSE, row.names = FALSE)
