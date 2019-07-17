

setwd("C:/Github/AIL_B6xBFMI/Candidate Genes")

myvcf <- read.csv("all_combined.vcf", sep = "\t", skip = 169, header=FALSE, colClasses="character")
dim(myvcf)

colnames(myvcf) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "BFMI860-12", "BFMI861-S1", "BFMI861-S2")

GTs1 <- unlist(lapply(strsplit(myvcf[,"BFMI861-S1"], ":"), "[", 1))
GTs1[which(GTs1 == "./.")] <- NA

interesting <- which(GTs1 != "0/0")
myvcf <- myvcf[interesting,]

complexSNPs <- grep(",", myvcf[,"ALT"])
myvcf <- myvcf[-complexSNPs,]

# Add a vep like 'location' column to the vcf file
myvcf <- cbind(location = as.character(paste0(myvcf[,"CHROM"], ":", myvcf[, "POS"], "-", myvcf[, "POS"])), myvcf)

#Keep only the non-duplicated entries, (remove the duplicated entries)
myvcf <- myvcf[which(!duplicated(myvcf[,"location"])),]

# Load the VEP results
myvep <- read.csv("vep_predictions_6_18_2019.txt", sep = "\t", skip = 1, header=FALSE, colClasses="character")
colnames(myvep) <- c("Uploaded_variation", "Location", "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID", "TSL", "APPRIS", "REFSEQ_MATCH", "SIFT", "CLIN_SIG", "SOMATIC", "PHENO", "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE")

allgenes <- unique(myvep[which(myvep[,"SYMBOL"] == "allgenes" & myvep[, "Consequence"] == "missense_variant"),"Location"])
myvcf[which(myvcf[, "location"] %in% allgenesmis),]

# Keep only the VEP entries that are in the VCF file
myvep <- myvep[which(myvep[, "Location"] %in% myvcf[, "location"]),]

# Do the genes first since we're looking for non-syn SNPs
vepgenes <- myvep[which(myvep[, "SYMBOL"] != "-"),]
allgenes <- unique(vepgenes[, "SYMBOL"])

resM <- matrix(NA, length(allgenes), 4)
resM[,1] <- allgenes
colnames(resM) <- c("name", "chr", "missense variant", "consequence")

mrow <- 1
for (gene in allgenes) {
  snpingene <- vepgenes[which(vepgenes[,"SYMBOL"] == gene),]
  chr <- unique(unlist(lapply(strsplit(snpingene[, "Location"], ":"), "[",1)))
  resM[mrow, "chr"] <- chr
  hasMissense <- any(unique(snpingene[,"Consequence"]) == "missense_variant")
  if(hasMissense) {
    resM[mrow, "missense"] <- "+"
    if(length(grep("tolerated", snpingene[,"SIFT"])) > 0) resM[mrow, "consequence"] <- "-"
    if(length(grep("tolerated_low_confidence", snpingene[,"SIFT"])) > 0) resM[mrow, "consequence"] <- "0"
    if(length(grep("deleterious_low_confidence", snpingene[,"SIFT"])) > 0) resM[mrow, "consequence"] <- "+"
    if(length(grep("deleterious", snpingene[,"SIFT"])) > 0) resM[mrow, "consequence"] <- "++"
  }else{
    resM[mrow, "missense"] <- "-"
  }
  mrow <- mrow + 1
}
resM[which(resM[,"consequence"] == "++"),]





for (gene in allgenes) {
  snpingene <- vepgenes[which(vepgenes[,"Expression Level"] == gene),]
  Expressioninoneofthe4regions <- any(unique(snpingene[,"Expression"]) == "brain", "liver", "skeleton", "hypothalamus")
  if(Expressioninoneofthe4regions)
  }













