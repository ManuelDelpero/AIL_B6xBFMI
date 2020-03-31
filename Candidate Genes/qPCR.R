# RT pcr analysis
#
# copyright (c) 2018-2021 - Brockmann group - HU Berlin, Danny Arends and Manuel Delpero
# first written february, 2020


setwd("C:/Users/Manuel/Desktop/AIL_B6xBFMI/RAWDATA")
mdata <- read.csv("qPCR_Sandra_Feb.txt", sep = "\t", na.strings = c("Undetermined", ""), colClasses = "character")
# Remove things with missing CT values
mdata <- mdata[-which(is.na(mdata[, "CT"])),]

# Fix some issues
mdata[,"Target.Name"] <- tolower(mdata[,"Target.Name"])
mdata[which(mdata[, "Sample.Name"] %in% c("primt0,1", "primt1", "primt10")),"Sample.Name"] <- "primertest"

# Our marker / Genotype data
topmarkerChr3 <- "gUNC5036315"
topmarkerChr8 <- "S1H083826428"

marker.data <- read.table("genomatrix.clean.txt", sep = "\t", colClasses = "character")
genotypes <- marker.data[c(topmarkerChr3, topmarkerChr8),]
colnames(genotypes) <- gsub("X666.", "", colnames(genotypes), fixed=TRUE)
genotypes <- genotypes[which(colnames(genotypes) %in% mdata[, "Sample.Name"])]

### Group into genotypes by the marker near the CES genes
mdata <- mdata[which(mdata[, "Sample.Name"] %in% colnames(genotypes)),]

mdata <- cbind(mdata, t(genotypes[,mdata[, "Sample.Name"]]))
B6_B6 <- which(mdata[, topmarkerChr3] != "AA" & mdata[, topmarkerChr8] != "AA")
BFMI_B6 <- which(mdata[, topmarkerChr3] == "AA" & mdata[, topmarkerChr8] != "AA")
X_BFMI <- which(mdata[, topmarkerChr8] == "AA")

mdata <- cbind(mdata, GT = NA)
mdata[B6_B6, "GT"] <- "B6_B6"
mdata[BFMI_B6, "GT"] <- "BFMI_B6"
mdata[X_BFMI, "GT"] <- "X_BFMI"

genotypes <- mdata[, "GT"]
names(genotypes) <- mdata[,"Sample.Name"]

samples <- unique(mdata[, "Sample.Name"])
genes <- unique(mdata[, "Target.Name"])
housekeeper <- c("actb", "rsp")
genes <- genes[which(!(genes %in% housekeeper))]

### Copy the original data, so we can compare back
mdataOld <- mdata

redo <- c() # samples that need to be redone
# Fix/Remove the large CT.SD values
for (g in c(genes,housekeeper)) {
  for (s in samples) {
    idx <- which(mdata[, "Sample.Name"] == s &  mdata[, "Target.Name"] == g)
    gData <- mdata[idx,]
    vals <- as.numeric(gData[, "CT"])
    valsO <- as.numeric(gData[, "CT"])
    CTmean <- round(median(vals, na.rm = TRUE), 1)
    CTsd <- round(sd(vals, na.rm = TRUE), 1)
    #if(is.na(CTsd)) stop("!!!")
    if (is.na(CTsd) || CTsd > 0.2) {
      cat(s, " ", g, " CTSD is too large: ", CTsd, ":", as.numeric(gData[, "CT"]), ", ")
      # If we have 3 values figure out if removing the one furthest from the mean makes the data suitable
      if (length(vals) >= 3) {
        CTsdC <- CTsd
        biggestDifs <- c()
        while(length(vals) > 2 && CTsdC > 0.2){
          biggestDif <- which.max(abs(vals - CTmean))
          vals <- vals[-biggestDif]
          CTsdC <- round(sd(vals, na.rm = TRUE), 1)
          biggestDifs <- c(biggestDifs, biggestDif)
        }
        if (CTsdC <= 0.2) {
          cat("Fix ",length(biggestDifs),", removed:",valsO[biggestDifs],"\n")
          mdata[idx[biggestDifs],"CT"] <- NA
        } else {
          cat("Removing ALL\n")
          redo <- rbind(redo, c(s,g))
          mdata[idx,"CT"] <- NA
        }
      } else {
        cat("No 3 measurements, Removing ALL\n")
        redo <- rbind(redo, c(s,g))
        mdata[idx,"CT"] <- NA
      }
    }
  }
}

write.table(redo, "failed_qPCR_toREDO.txt", sep = "\t", quote = FALSE, row.names=FALSE)

# Compute CT values relative to the housekeepers
pvals <- c()
for(g in genes){
  mymatrix <- c()
  for(s in samples){
    gData <- mdata[which(mdata[, "Sample.Name"] == s &  mdata[, "Target.Name"] == g),]
    gMean <- round(mean(as.numeric(gData[, "CT"]), na.rm=TRUE),1)
    myrow <- c(g, s, genotypes[s], gMean)
    for(hk in housekeeper){
      hKeeper <- mdata[which(mdata[, "Sample.Name"] == s &  mdata[, "Target.Name"] == hk),]
      hKeeperMean <- round(mean(as.numeric(hKeeper[, "CT"]), na.rm=TRUE),1)
      myrow <- c(myrow, hKeeperMean)
    }
    mymatrix <- rbind(mymatrix, myrow)
  }
  rownames(mymatrix) <- samples
  colnames(mymatrix) <- c("Gene", "Sample", "GT", "CT", housekeeper)
  dCTlist <- vector("list", length(unique(mymatrix[, "GT"])))
  names(dCTlist) <- unique(mymatrix[, "GT"])
  for(gt in unique(mymatrix[, "GT"])){
    gtmatrix <- mymatrix[which(mymatrix[, "GT"] == gt),]
    CTdiff_actb <- as.numeric(gtmatrix[, "CT"]) - as.numeric(gtmatrix[, "actb"])
    CTdiff_actb <- 2^(-CTdiff_actb)
    CTdiff_rsp <- as.numeric(gtmatrix[, "CT"]) - as.numeric(gtmatrix[, "rsp"])
    CTdiff_rsp <- 2^(-CTdiff_rsp)
    mm <- cbind(CTdiff_actb, CTdiff_rsp)
    rownames(mm) <- rownames(gtmatrix)
    dCTlist[[gt]] <- mm
  }
  # Adjust the values setting the B6 group to 1
  act <- mean(dCTlist$B6_B6[,1], na.rm = TRUE)
  rsp <- mean(dCTlist$B6_B6[,2], na.rm = TRUE)
  dCTlist$B6_B6[,1] <- dCTlist$B6_B6[,1] / act
  dCTlist$B6_B6[,2] <- dCTlist$B6_B6[,2] / rsp
  dCTlist$X_BFMI[,1] <- dCTlist$X_BFMI[,1] / act
  dCTlist$X_BFMI[,2] <- dCTlist$X_BFMI[,2] / rsp
  dCTlist$BFMI_B6[,1] <- dCTlist$BFMI_B6[,1] / act
  dCTlist$BFMI_B6[,2] <- dCTlist$BFMI_B6[,2] / rsp
  for(gtA in c("X_BFMI", "B6_B6")){
    for(gtB in c("B6_B6", "BFMI_B6")){
      if(gtA != gtB){
        tryCatch({
        P_actb <- t.test(dCTlist[[gtA]][, "CTdiff_actb"], dCTlist[[gtB]][, "CTdiff_actb"])$p.value
        gtAgroup_actb <- paste0(names(which(!is.na(dCTlist[[gtA]][, "CTdiff_actb"]))), collapse=", ")
        gtBgroup_actb <- paste0(names(which(!is.na(dCTlist[[gtB]][, "CTdiff_actb"]))), collapse=", ")
        
        P_rsp <- t.test(dCTlist[[gtA]][, "CTdiff_rsp"], dCTlist[[gtB]][, "CTdiff_rsp"])$p.value
        gtAgroup_rsp <- paste0(names(which(!is.na(dCTlist[[gtA]][, "CTdiff_rsp"]))), collapse=", ")
        gtBgroup_rsp <- paste0(names(which(!is.na(dCTlist[[gtB]][, "CTdiff_rsp"]))), collapse=", ")
        
        pvals <- rbind(pvals, c(g, gtA, gtB, "actb", round(mean(dCTlist[[gtA]][, "CTdiff_actb"],na.rm=TRUE),3), round(sd(dCTlist[[gtA]][, "CTdiff_actb"],na.rm=TRUE),3), round(mean(dCTlist[[gtB]][, "CTdiff_actb"],na.rm=TRUE),3), round(sd(dCTlist[[gtB]][, "CTdiff_actb"],na.rm=TRUE),3), P_actb, gtAgroup_actb, gtBgroup_actb))
        pvals <- rbind(pvals, c(g, gtA, gtB, "rsp", round(mean(dCTlist[[gtA]][, "CTdiff_rsp"],na.rm=TRUE),3), round(sd(dCTlist[[gtA]][, "CTdiff_rsp"],na.rm=TRUE),3), round(mean(dCTlist[[gtB]][, "CTdiff_rsp"],na.rm=TRUE),3), round(sd(dCTlist[[gtB]][, "CTdiff_rsp"],na.rm=TRUE),3), P_rsp, gtAgroup_rsp, gtBgroup_rsp))
        }, error = function(err) {
        
        })
      }
    }
  }
}
colnames(pvals) <- c("Gene", "GT1", "GT2", "Housekeeper", "Mean GT1", "SD GT1", "Mean GT2", "SD GT2", "p.value", "Ind GT1", "Ind GT2")
data.frame(pvals)

colnames(redo) <- c("Sample", "Gene")
redo

# Barplots to represent the results
genesCAP <- rbind(as.numeric(pvals[35, "Mean GT1"]), as.numeric(pvals[33, "Mean GT1"]), as.numeric(pvals[33, "Mean GT2"])) 
genesCAP <- cbind(genesCAP, c("Capns2[B6_B6]", "Capns2[x_BFMI]", "Capns2[BFMI_B6]"))
genes1B <- rbind(as.numeric(pvals[2,"Mean GT2"]), as.numeric(pvals[2, "Mean GT1"]), as.numeric(pvals[4, "Mean GT2"])) 
genes1B <- cbind(genes1B, c("Ces1b[B6_B6]", "Ces1b[x_BFMI]", "Ces1b[BFMI_B6]"))
genes1C <- rbind(as.numeric(pvals[8, "Mean GT2"]), as.numeric(pvals[8, "Mean GT1"]), as.numeric(pvals[10, "Mean GT2"])) 
genes1C <- cbind(genes1C, c("Ces1b[B6_B6]", "Ces1b[x_BFMI]", "Ces1b[BFMI_B6]"))
genes1D <- rbind(as.numeric(pvals[13, "Mean GT2"]), as.numeric(pvals[13, "Mean GT1"]), as.numeric(pvals[15, "Mean GT2"])) 
genes1D <- cbind(genes1D, c("Ces1D[B6_B6]", "Ces1D[x_BFMI]", "Ces1D[BFMI_B6]"))
genes1E <- rbind(as.numeric(pvals[24, "Mean GT1"]), as.numeric(pvals[22, "Mean GT1"]), as.numeric(pvals[24, "Mean GT2"])) 
genes1E <- cbind(genes1E, c("Ces1E[B6_B6]", "Ces1E[x_BFMI]", "Ces1E[BFMI_B6]"))
genes1F <- rbind(as.numeric(pvals[25, "Mean GT2"]), as.numeric(pvals[25, "Mean GT1"]), as.numeric(pvals[29, "Mean GT2"])) 
genes1F <- cbind(genes1F, c("Ces1F[B6_B6]", "Ces1F[x_BFMI]", "Ces1F[BFMI_B6]"))
genes <- cbind(genesCAP, genes1B, genes1C, genes1D, genes1E, genes1F)

STCAP <- rbind(as.numeric(pvals[35, "SD GT1"]), as.numeric(pvals[33, "SD GT1"]), as.numeric(pvals[33, "SD GT2"])) 
STCAP <- cbind(STCAP, c("Capns2[B6_B6]", "Capns2[x_BFMI]", "Capns2[BFMI_B6]"))
ST1B <- rbind(as.numeric(pvals[2,"SD GT2"]), as.numeric(pvals[2, "SD GT1"]), as.numeric(pvals[4, "SD GT2"])) 
ST1B <- cbind(ST1B, c("Ces1b[B6_B6]", "Ces1b[x_BFMI]", "Ces1b[BFMI_B6]"))
ST1C <- rbind(as.numeric(pvals[8, "SD GT2"]), as.numeric(pvals[8, "SD GT1"]), as.numeric(pvals[10, "SD GT2"])) 
ST1C <- cbind(ST1C, c("Ces1b[B6_B6]", "Ces1b[x_BFMI]", "Ces1b[BFMI_B6]"))
ST1D <- rbind(as.numeric(pvals[13, "SD GT2"]), as.numeric(pvals[13, "SD GT1"]), as.numeric(pvals[15, "SD GT2"])) 
ST1D <- cbind(ST1D, c("Ces1D[B6_B6]", "Ces1D[x_BFMI]", "Ces1D[BFMI_B6]"))
ST1E <- rbind(as.numeric(pvals[24, "SD GT1"]), as.numeric(pvals[22, "SD GT1"]), as.numeric(pvals[24, "SD GT2"])) 
ST1E <- cbind(ST1E, c("Ces1E[B6_B6]", "Ces1E[x_BFMI]", "Ces1E[BFMI_B6]"))
ST1F <- rbind(as.numeric(pvals[25, "SD GT2"]), as.numeric(pvals[25, "SD GT1"]), as.numeric(pvals[29, "SD GT2"])) 
ST1F <- cbind(ST1F, c("Ces1F[B6_B6]", "Ces1F[x_BFMI]", "Ces1F[BFMI_B6]"))
SDgenes <- cbind(STCAP, ST1B, ST1C, ST1D, ST1E, ST1F)

colnames(SDgenes) <- c("Capns2", "Capns2geno", "Ces1b", "Ces1bgeno", "Ces1c", "Ces1cgeno", "Ces1d", "Ces1dgeno", "Ces1e", "Ces1egeno", "Ces1f", "Ces1fgeno")
colnames(genes) <- c("meanCapns2", "Capns2geno", "meanCes1b", "Ces1bgeno", "meanCes1c", "Ces1cgeno", "meanCes1d", "Ces1dgeno", "meanCes1e", "Ces1egeno", "meanCes1f", "Ces1fgeno")
mat <- cbind(as.numeric(genes[,1]), as.numeric(genes[,3]), as.numeric(genes[,5]), as.numeric(genes[,7]), as.numeric(genes[,9]), as.numeric(genes[,11]))

# Ces genes
x <- barplot(mat[,-1],beside=TRUE,col= c("gray20", "gray50", "gray88"), main = "Ces genes relative expression", ylab = "Relative gene expression", ylim = c(0, 8))
  axis(1, at = c(2.5, 6.5, 10.5, 14.5, 18.5), c("Ces1b", "Ces1c", "Ces1d", "Ces1e", "Ces1f"))
  y <- 1.7
  # set an offset for tick lengths
  offset <- 0.1
  # draw first horizontal line
  lines(x[2:3],c(y, y))
  lines(x[1:3],c(y + offset + offset, y + offset + offset, y + offset + offset))
  # draw asterics
  text(3,y+offset,"**")
  text(2.5,y+offset + offset + offset,"**")
  segments(x[,1], mat[,2] - as.numeric(SDgenes[,3]), x[,1], mat[,2] + as.numeric(SDgenes[,3]), lwd = 1.5)
  arrows(x[,1], mat[,2] - as.numeric(SDgenes[,3]), x[,1], mat[,2] + as.numeric(SDgenes[,3]), lwd = 1.5, angle = 90, code = 3, length = 0.05)
  segments(x[,2], mat[,3] - as.numeric(SDgenes[,5]), x[,2], mat[,3] + as.numeric(SDgenes[,5]), lwd = 1.5)
  arrows(x[,2], mat[,3] - as.numeric(SDgenes[,5]), x[,2], mat[,3] + as.numeric(SDgenes[,5]), lwd = 1.5, angle = 90, code = 3, length = 0.05)
  segments(x[,3], mat[,4] - as.numeric(SDgenes[,7]), x[,3], mat[,4] + as.numeric(SDgenes[,7]), lwd = 1.5)
  arrows(x[,3], mat[,4] - as.numeric(SDgenes[,7]), x[,3], mat[,4] + as.numeric(SDgenes[,7]), lwd = 1.5, angle = 90, code = 3, length = 0.05)
  segments(x[,4], mat[,5] - as.numeric(SDgenes[,9]), x[,4], mat[,5] + as.numeric(SDgenes[,9]), lwd = 1.5)
  arrows(x[,4], mat[,5] - as.numeric(SDgenes[,9]), x[,4], mat[,5] + as.numeric(SDgenes[,9]), lwd = 1.5, angle = 90, code = 3, length = 0.05)
  segments(x[,5], mat[,6] - as.numeric(SDgenes[,11]), x[,5], mat[,6] + as.numeric(SDgenes[,11]), lwd = 1.5)
  arrows(x[,5], mat[,6] - as.numeric(SDgenes[,11]), x[,5], mat[,6] + as.numeric(SDgenes[,11]), lwd = 1.5, angle = 90, code = 3, length = 0.05)
  legend("topright", legend=c("B6 Ch3/B6 Chr 8", "B6 or BFMI Chr3/BFMI Chr8", "BFMI Chr3/B6 Chr8"), fill= c("gray20", "gray50", "gray88"), bty = "n")
  
# Capns2
x <- barplot(mat[,1],beside=TRUE,col= c("gray20", "gray50", "gray88"), main = "Capns2 relative expression", ylab = "Relative gene expression", ylim = c(0, 1000))
  axis(1, at = 1.9, "Capns2")
  y <- 800
  # set an offset for tick lengths
  offset <- 20
  # draw first horizontal line
  lines(x[2:3],c(y, y))
  lines(x[1:3],c(y + offset + offset, y + offset + offset, y + offset + offset))
  # draw asterics
  text(2.5,y+offset,"**")
  text(1.9,y+offset + offset + offset,"**")
  segments(x, mat[,1] - as.numeric(SDgenes[,1]), x, mat[,1] + as.numeric(SDgenes[,1]), lwd = 1.5)
  arrows(x, mat[,1] - as.numeric(SDgenes[,1]), x, mat[,1] + as.numeric(SDgenes[,1]), lwd = 1.5, angle = 90, code = 3, length = 0.05)
  legend("topright", legend=c("B6 Ch3/B6 Chr 8", "B6 or BFMI Chr3/BFMI Chr8", "BFMI Chr3/B6 Chr8"), fill= c("gray20", "gray50", "gray88"), bty = "n")
