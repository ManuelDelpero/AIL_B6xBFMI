setwd("C:/Users/Manuel/Desktop/AIL_B6xBFMI/RAWDATA")

phenotypes <- read.csv("allPhenotypes.txt", header = TRUE, check.names = FALSE, sep="\t")
genotypes <- read.csv("genomatrix.clean.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")

topmarkerChr3 <- "gUNC5036315"
topmarkerChr8 <- "S1H083826428"

#AA = BFMI
markerdata <- cbind(topmarkerChr3 = unlist(genotypes[topmarkerChr3,]), 
                    topmarkerChr8 = unlist(genotypes[topmarkerChr8,]))

combineddata <- cbind(markerdata, TriGly = phenotypes[rownames(markerdata), "Triglycerides/Proteins"],
                                  d98    = phenotypes[rownames(markerdata), "d98"])

combineddata <- combineddata[-which(is.na(combineddata[,"TriGly"])),]

# BFMI on jObes1
BFMIlike <- combineddata[which(combineddata[, "topmarkerChr3"] == "AA"),]
BFMI_BFMI <- BFMIlike[which(BFMIlike[, "topmarkerChr8"] == "AA"),]
BFMI_BFMI <- BFMI_BFMI[sort(as.numeric(BFMI_BFMI[, "TriGly"]), index.return=TRUE, decreasing=TRUE)$ix,]

BFMI_B6 <- BFMIlike[which(BFMIlike[, "topmarkerChr8"] != "AA"),]
BFMI_B6 <- BFMI_B6[sort(as.numeric(BFMI_B6[, "TriGly"]), index.return=TRUE)$ix,]

# B6 or HET on jObes1
B6like <- combineddata[which(combineddata[, "topmarkerChr3"] != "AA"),]

B6_BFMI <- B6like[which(B6like[, "topmarkerChr8"] == "AA"),]
B6_BFMI <- B6_BFMI[sort(as.numeric(B6_BFMI[, "TriGly"]), index.return=TRUE, decreasing=TRUE)$ix,]

B6_B6 <- B6like[which(B6like[, "topmarkerChr8"] != "AA"),]
B6_B6 <- B6_B6[sort(as.numeric(B6_B6[, "TriGly"]), index.return=TRUE)$ix,]

# Do the t-test to compare the groups
t.test(as.numeric(B6_B6[, "TriGly"]), as.numeric(BFMI_B6[, "TriGly"]))
t.test(as.numeric(B6_BFMI[, "TriGly"]), as.numeric(BFMI_BFMI[, "TriGly"]))

#Result, chromosome 3 matters if you're a B6 on chromosome 8
# if you're BFMI on chromosome 8 then your chromosome 3 status doesn't matters
# As such we can define 3 groups: B6_B6, BFMI_B6, (B6_BFMI and BFMI_BFMI)
# We select some animals from each of them !

BFMI_C8 <- rbind(B6_BFMI, BFMI_BFMI)
BFMI_C8 <- BFMI_C8[sort(as.numeric(BFMI_C8[, "TriGly"]), index.return=TRUE, decreasing=TRUE)$ix,]

#Define our groups (select 20 in total)
G1 <- rownames(B6_B6[1:6,]) #6 from B6 B6
G2 <- rownames(BFMI_B6[1:6,]) # 6 from BFMI B6
G3 <- rownames(BFMI_C8[2:9,]) # 8 from BFMI_C8, from 2 since #1 is an outlier

days <- colnames(phenotypes)[7:24]
numdays <- as.numeric(gsub("d", "", days))

pdf("GrowthCurves_Selected")
### Create the plot
plot(x = c(min(numdays), max(numdays)), y = c(10, 55), t = 'n', ylab="Bodyweight [g]", xlab="Time [days]")
for(i in G1){
  bfmiS <- as.numeric(combineddata[i, "topmarkerChr3"] == "AA")
  points(numdays, phenotypes[i, days], col = "red", t = 'b', pch = 16 + bfmiS)
}
for(i in G2){
  bfmiS <- as.numeric(combineddata[i, "topmarkerChr3"] == "AA")
  points(numdays, phenotypes[i, days], col = "blue", t = 'b', pch = 16 + bfmiS)
}
for(i in G3){
  bfmiS <- as.numeric(combineddata[i, "topmarkerChr3"] == "AA")
  points(numdays, phenotypes[i, days], col = "orange", t = 'b', pch = 16 + bfmiS)
}
legend("topleft", c("B6 B6", "BFMI B6", "??? BFMI", "jObes1[b6]", "jObes1[BFMI]"), 
                  col=c("red", "blue", "orange", "black", "black"), lwd=c(1,1,1,NA,NA), pch =c(NA,NA,NA, 16, 17))
				
dev.off()
