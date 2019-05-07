#
# Use biomaRt to download all the protein coding genes in a several QTL regions
#

# If not installed, install biomart by using:
#   install.packages("BiocManager")
#   BiocManager::install("biomaRt")

library(biomaRt)

setwd("D:/Ddrive/Collegues/Sandra")

regions <- read.table("input.matrix.txt", sep = "\t", header = TRUE)

bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

for (mrow in 1:nrow(regions)) {
  region <- paste0(regions[mrow,"Chromosome"], ":", regions[mrow,"Start"], ":", regions[mrow,"End"])
  cat("line: ", mrow, " has region: ", region, "\n")
  res.biomart <- getBM(attributes = c("ensembl_gene_id",                                                    # Things that we want to get from biomart
                                      "chromosome_name", "start_position", "end_position", "strand", 
                                      "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                       filters = c("chromosomal_region", "biotype"),                                        # Things that we will use to query biomart
                       values = list(region, "protein_coding"),                                             # The thing that we are querying
                       mart = bio.mart)
  filename <- paste0(regions[mrow,"Phenotype"], "_", regions[mrow,"Chromosome"], "-", regions[mrow,"Start"], "_", regions[mrow,"End"], ".genes")
  cat("line: ", mrow, " has ", nrow(res.biomart), "genes, will save to file:", filename, "\n")
  write.table(res.biomart, file=filename, sep="\t")
}
