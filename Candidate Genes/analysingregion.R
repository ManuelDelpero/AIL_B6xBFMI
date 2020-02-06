#written by Manuel Delpero
#first written May, 2019

setwd("/home/manuel/AIL_B6xBFMI/RAWDATA_GENERAL/SNPs/")
regions <- read.table("QTL_regions_Lod1.5.txt", sep = "\t", header = TRUE)

# Get genes in regions
genes <- vector("list", nrow(regions))
for(x in 1:nrow(regions)){																				#for loop number of rows in the region
  # check if the region is valid (not NA)
  if(!is.na(regions[x, "Chr"])){
    # get the genes in the region
    genes[[x]] <- getregion(bio.mart, regions[x, "Chr"], regions[x, "StartPos"], regions[x, "StopPos"])
	cat(x, " has ", nrow(genes[[x]]), "genes\n")
    fname <- paste0("genes_in_", regions[x, "Chr"],"-", regions[x, "StartPos"], ":", regions[x, "StopPos"], ".txt")
    write.table(genes[[x]], file = fname, sep="\t", quote = FALSE)
  }else{
 	cat(x, " has NA region\n")
  }
}

# figure out all the unique genes, result: matrix with 4 columns: name, chromosome, start position, end position
uniquegenes <- NULL
for(x in genes){ 
  if(!is.null(x)){
    subset <- x[ ,c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand")]
	uniquegenes <- rbind(uniquegenes, subset) 
  }
}
uniquegenes <- uniquegenes[!duplicated(uniquegenes),] 
table(uniquegenes[ ,"chromosome_name"])

bamfiles <- c("/halde/BFMI_Alignment_Mar19/merged_sorted_860-S12.bam",  # 860-S12  (high coverage)
             "/halde/BFMI_Alignment_Mar19/merged_sorted_861-S1.bam",    # 861-S1 (medium coverage)
             "/halde/BFMI_Alignment_Mar19/merged_sorted_861-S2.bam")    # 861-S2 (medium coverage)

# Snps in genes
for(x in 1:nrow(uniquegenes)){ 
  startpos <- uniquegenes[x, 3]
  endpos <- uniquegenes[x, 4]
  if(uniquegenes[x, 5] == 1){
    startpos <- startpos-500
  }else{
    endpos <- endpos +500
  }
  callSNPs(bamfiles, uniquegenes[x, 2], startpos, endpos, uniquegenes[x, 1]) 
}


filelist <- list.files(".") #we removed all the .txt files

allSNPs <- NULL
for(file in filelist){
  if(length(readLines(file)) > 169){
    mcontent <- read.csv(file, sep = "\t", skip = 169, header=FALSE, colClasses=c("character"))
    gene = gsub(".snps-filtered.vcf", "", file, fixed=T)
    allSNPs <- rbind(allSNPs, cbind(gene, mcontent))
  }
}

# Sort by chromosomes and position otherwise VEP is complaining
allSNPs <- allSNPs[order(as.numeric(allSNPs[,3])),]
chromosomes <- c(1, 3, 5, 8)

annotation <- c()
for(chr in chromosomes){
  annotation <- rbind(annotation, allSNPs[as.numeric(allSNPs[,"V1"]) == chr,])
}

allSNPs <- annotation

header = readLines(filelist[1], n = 169)
cat(paste0(header, collapse = "\n"), "\n", file = "all_combined.vcf")
# File containing all SNPs in the genes
write.table(allSNPs[,-1], file = "all_combined.vcf", sep = "\t", quote=FALSE, append = TRUE, col.names=FALSE, row.names= FALSE)