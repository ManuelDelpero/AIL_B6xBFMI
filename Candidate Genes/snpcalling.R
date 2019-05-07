#
# Use bcftools to call SNPs from bam files
#

execute <- function(x, intern = FALSE, execute = TRUE){
  cat("----", x, "\n")
  if (execute) {
    res <- system(x, intern = intern)
    cat(">>>>", res[1], "\n")
    if(res[1] >= 1) q("no")
  }
}

callSNPs <- function(bamfiles, chr = 1, startpos = 1, endpos = 2, outname = "mySNPs") {
  bamstr = paste0(bamfiles, collapse = " ") # Collapse all the bam files in a single string
  bcftools = "/home/danny/Github/bcftools/bcftools" # Location of BCFtools executable
  reference = "/home/danny/References/Mouse/GRCm38_95/Mus_musculus.GRCm38.dna.toplevel.fa" #Reference genome
  region = paste0(chr, ":", format(startpos, scientific = FALSE), "-", format(endpos, scientific = FALSE)) # Region requested

  cmd1 <- paste0(bcftools, " mpileup -q 30 -Ou -r ", region, " -f ", reference, " ", bamstr)
  cmd2 <- paste0(bcftools, " call -mv -Ov ")
  cmd3 <- paste0(bcftools, " view -v snps -i '%QUAL>=20 && DP>10' - -o ", outname, ".snps-filtered.vcf")
  execute(paste0(cmd1, " | ", cmd2, " | ", cmd3))
  invisible("")
}

callSNPs(c("/halde/BFMI_Alignment_Mar19/merged_chr19.bam", #860-S12 high coverage
           "/home/sandra/NAS/DNA/Sequencing/DifeMouse/RAW/BAM/ext_L7256.bam", # 861-S1 (Needs realignment)
           "/home/sandra/NAS/DNA/Sequencing/DifeMouse/RAW/BAM/ext_L7257.bam" # 861-S2 (Needs realignment)
           ), 19, 3000000, 3500000, "test")
