# AIL_B6xS1 
#
# copyright (c) - Danny Arends & Manuel Delpero
# first written september, 2019
# LMMMQM for bodyweight using AIL_B6xS1

setwd("C:/Users/Manuel/Desktop/AIL_B6xBFMI/RAWDATA")

genotypes <- read.csv("genomatrix.clean_numeric.txt", header = TRUE, check.names = FALSE, sep="\t", colClasses="character")
phenotypes <- read.csv("allPhenotypes.txt", header = TRUE, check.names = FALSE, sep="\t", row.names=1)
markerannot <- read.csv("SNP_Map.txt", header=TRUE, sep="\t", row.names=2, check.names=FALSE)
markerannot <- markerannot[,-1]
markerannot <- markerannot[rownames(genotypes),]
markerannot <- markerannot[sort(markerannot[,"Position"], index.return=TRUE)$ix,]

chromosomes <- c(1:19, "X", "Y")

annotation <- c()
for(chr in chromosomes){
  annotation <- rbind(annotation, markerannot[markerannot[,"Chromosome"] == chr,])
}

phenonames <- colnames(phenotypes)[-c(1:5)]

# Make sure that the ordering between phenotypes and genotypes matches !!!!!
# Also sort the markers by their natural chromosome ordering
genotypes <- genotypes[rownames(annotation), rownames(phenotypes)]
write.table(genotypes, "OrderedGenotypes.txt", sep = "\t", quote=FALSE)

library(lme4)
ctrl = lmerControl(optimizer ="Nelder_Mead", optCtrl=list(maxfun=20000))

# Covariates we need to include in the model / tested previously
individual <- c()
timepoint <- c()
weight <- c()
sex <- c()
wg <- c()
mother <- c()
for(tp in paste0("d", seq(21, 140, 7))){
  individual <- c(individual, rownames(phenotypes))
  timepoint <- c(timepoint, as.numeric(rep(gsub("d", "", tp), nrow(phenotypes))))
  weight <- c(weight, phenotypes[, tp])
  sex <- c(sex, as.character(phenotypes[, "Sex"]))
  wg <- c(wg, phenotypes[, "WG"])
  mother <- c(mother, as.character(phenotypes[, "Mutter"]))
}

### Forward modeling figure out what needs to be in the model
model1 <- lmer(weight ~ (1|individual), REML=FALSE, control = ctrl)
model2 <- lmer(weight ~ timepoint + (1|individual), REML=FALSE, control = ctrl)
AIC(model1, model2)
model2b <- lmer(weight ~ timepoint + (timepoint|individual), REML=FALSE, control = ctrl)
AIC(model2, model2b)

model4 <- lmer(weight ~ sex + timepoint + (timepoint|individual), REML=FALSE, control = ctrl)
AIC(model3, model4)
model5 <- lmer(weight ~ mother + sex + timepoint + (timepoint|individual), REML=FALSE, control = ctrl)
AIC(model4, model5)
model6 <- lmer(weight ~ mother + sex + timepoint + I(timepoint^2) + (timepoint|individual), REML=FALSE, control = ctrl)
AIC(model5, model6)
model7 <- lmer(weight ~ mother + sex + timepoint + I(timepoint^2) + I(timepoint^3) + (timepoint|individual), REML=FALSE, control = ctrl)
AIC(model6, model7)


mpm <- cbind(individual, timepoint, weight, sex, wg, mother)

### Mapping across the genome using the extended model found above
# TODO put in chromosome 3 and 6, drop 3 when close to Bbs7, drop 6 when close to topmarker chr 6
resM <- c()
for(x in 1:nrow(annotation)){
  mgt.AD <- as.numeric(unlist(genotypes[x,mpm[,1]]))
  mgt.DD <- as.numeric(as.numeric(unlist(genotypes[x,mpm[,1]])) != 0)
  mgt.topChr6 <- as.numeric(unlist(genotypes["gUNC10595065",mpm[,1]]))
  mgt.topChr3 <- as.numeric(unlist(genotypes["gUNC5046545", mpm[,1]]))
  
  isMissing <- c(which(is.na(mgt.AD)), which(is.na(mgt.DD)))

  weight <- as.numeric(mpm[,"weight"])
  individual <- mpm[,"individual"]
  timepoint <- as.numeric(mpm[,"timepoint"])
  sex <- mpm[,"sex"]
  wg <- as.numeric(mpm[,"wg"])
  mother <- mpm[,"mother"]
  if(length(isMissing) > 0){
    weight <- weight[-isMissing]
    individual <- individual[-isMissing]
    timepoint <- timepoint[-isMissing]
    sex <- sex[-isMissing]
    wg <- wg[-isMissing]
    mother <- mother[-isMissing]
    mgt.AD <- mgt.AD[-isMissing]
    mgt.DD <- mgt.DD[-isMissing]
  }
  
  if ((annotation[x, "Chromosome"] == 3) && (annotation[x, "Position"] >=  annotation["gUNC5046545","Position"] - 5000000) && (annotation[x, "Position"] <= annotation["gUNC5046545","Position"] + 5000000)){
    tryCatch(
      model.full <- lmer(weight ~ mother + sex + timepoint + I(timepoint^2) + I(timepoint^3) + (timepoint|individual) + mgt.topChr6 + mgt.AD + mgt.AD:timepoint + mgt.DD + mgt.DD:timepoint , REML=FALSE, control = ctrl)
      , error = function(e) e)
    tryCatch(
      model.add <- lmer(weight ~ mother + sex + timepoint + I(timepoint^2) + I(timepoint^3) + (timepoint|individual) + mgt.topChr6 + mgt.AD + mgt.AD:timepoint, REML=FALSE, control = ctrl)
      , error = function(e) e)
    tryCatch(
      model.null <- lmer(weight ~ mother + sex + timepoint + I(timepoint^2) + I(timepoint^3) + (timepoint|individual) + mgt.topChr6, REML=FALSE, control = ctrl)
      , error = function(e) e)
  } else if ((annotation[x, "Chromosome"] == 6) && (annotation[x, "Position"] >=  annotation["gUNC10595065","Position"] - 5000000) && (annotation[x, "Position"] <= annotation["gUNC10595065","Position"] + 5000000)){
    tryCatch(
      model.full <- lmer(weight ~ mother + sex + timepoint + I(timepoint^2) + I(timepoint^3) + (timepoint|individual) + mgt.topChr3 + mgt.AD + mgt.AD:timepoint + mgt.DD + mgt.DD:timepoint, REML=FALSE, control = ctrl)
      , error = function(e) e)
    tryCatch(
      model.add <- lmer(weight ~ mother + sex + timepoint + I(timepoint^2) + I(timepoint^3) + (timepoint|individual) + mgt.topChr3 + mgt.AD + mgt.AD:timepoint, REML=FALSE, control = ctrl)
      , error = function(e) e)
    tryCatch(
      model.null <- lmer(weight ~ mother + sex + timepoint + I(timepoint^2) + I(timepoint^3) + (timepoint|individual) + mgt.topChr3, REML=FALSE, control = ctrl)
      , error = function(e) e)
  } else{
   tryCatch(
     model.full <- lmer(weight ~ mother + sex + timepoint + I(timepoint^2) + I(timepoint^3) + (timepoint|individual) + mgt.topChr3 + mgt.topChr6 + mgt.AD + mgt.AD:timepoint + mgt.DD + mgt.DD:timepoint + mgt.topChr3 + mgt.topChr6, REML=FALSE, control = ctrl)
     , error = function(e) e)
   tryCatch(
     model.add <- lmer(weight ~ mother + sex + timepoint + I(timepoint^2) + I(timepoint^3) + (timepoint|individual) + mgt.topChr3 + mgt.topChr6 + mgt.AD + mgt.AD:timepoint, REML=FALSE, control = ctrl)
     , error = function(e) e)
   tryCatch(
     model.null <- lmer(weight ~ mother + sex + timepoint + I(timepoint^2) + I(timepoint^3) + (timepoint|individual) + mgt.topChr3 + mgt.topChr6, REML=FALSE, control = ctrl)
     , error = function(e) e)
  }
	 
  pval <- rep(NA, 3)
  tryCatch(pval[1] <- anova(model.full, model.add)["model.full", "Pr(>Chisq)"], error = function(e) e) # Dominance
  tryCatch(pval[2] <- anova(model.add, model.null)["model.add", "Pr(>Chisq)"], error = function(e) e) # Additive
  tryCatch(pval[3] <- anova(model.full, model.null)["model.full", "Pr(>Chisq)"], error = function(e) e) # Dom + Add
  resM <- rbind(resM, pval)
}

write.table(resM, file = "PvalLMMMQM.txt", sep = "\t", row.names = FALSE)
annotresM <- cbind(annotation, resM)

# Lod curve for Chr3
plot(x = c(1,nrow(annotresM[which(annotresM[,1] == "3"),])), y = c(0, 16))
lines(-log10(annotresM[which(annotresM[,1] == "3"), "3"]), col = "blue")
lines(-log10(annotresM[which(annotresM[,1] == "3"), "1"]), col = "black")
lines(-log10(annotresM[which(annotresM[,1] == "3"), "2"]), col = "orange")
  abline(h = -log10(0.05 / nrow(resM)), col="red")
  abline(h = -log10(0.01 / nrow(resM)), col="green")
  
# Manhattan plot
plot(x = c(1,nrow(annotresM)), y = c(0, 16))
lines(-log10(annotresM[, "1"]), col = "blue")
lines(-log10(annotresM[, "2"]), col = "black")
lines(-log10(annotresM[, "3"]), col = "orange")
  abline(h = -log10(0.05 / nrow(resM)), col="red")
  abline(h = -log10(0.01 / nrow(resM)), col="green")

topmarker <- "gUNC5046545"
mgt.AD <- as.numeric(unlist(genotypes[topmarker, mpm[,1]]))
mgt.DD <- as.numeric(as.numeric(unlist(genotypes[topmarker, mpm[,1]])) != 0)
isMissing <- c(which(is.na(mgt.AD)), which(is.na(mgt.DD)))

weight <- as.numeric(mpm[,"weight"])
individual <- mpm[,"individual"]
timepoint <- as.numeric(mpm[,"timepoint"])
sex <- mpm[,"sex"]
wg <- as.numeric(mpm[,"wg"])
mother <- mpm[,"mother"]
if(length(isMissing) > 0){
  weight <- weight[-isMissing]
  individual <- individual[-isMissing]
  timepoint <- timepoint[-isMissing]
  sex <- sex[-isMissing]
  wg <- wg[-isMissing]
  mother <- mother[-isMissing]
  mgt.AD <- mgt.AD[-isMissing]
  mgt.DD <- mgt.DD[-isMissing]
}

lmemodel <- lmer(weight ~ mother + sex + timepoint + I(timepoint^2) + I(timepoint^3) + (timepoint|individual) + mgt.AD + mgt.AD:timepoint + mgt.DD + mgt.DD:timepoint, REML=FALSE, control = ctrl)
fixedeffects <- fixef(lmemodel)
randomeffects <- ranef(lmemodel)$individual
effectmatrix <- matrix(NA, length(rownames(phenotypes)), 12, dimnames=list(rownames(phenotypes), c("fIntercept", "Mother", "Sex", "tp", "tp2", "tp3", "rIntercept", "rTp", "marker", "markertime", "markerDD", "markertimeDD")))
for(x in rownames(phenotypes)){
    eff.fixed.intercept <- fixedeffects["(Intercept)"]
    effectmatrix[x, 1] <- eff.fixed.intercept

    ind.mother <- as.character(phenotypes[x, "Mutter"])
    eff.mother <- fixedeffects[paste0("mother", ind.mother)]
    if (is.na(eff.mother)) eff.mother <- 0
    effectmatrix[x, 2] <- eff.mother
    
    ind.sex <- as.character(phenotypes[x, "Sex"])
    eff.sex <- fixedeffects[paste0("sex", ind.sex)]
    if (is.na(eff.sex)) eff.sex <- 0
    effectmatrix[x, 3] <- eff.sex
    
    eff.fixed.timepoint <- fixedeffects["timepoint"]
    effectmatrix[x, 4] <- eff.fixed.timepoint
    
    eff.fixed.timepoint2 <- fixedeffects["I(timepoint^2)"]
    effectmatrix[x, 5] <- eff.fixed.timepoint2
    
    eff.fixed.timepoint3 <- fixedeffects["I(timepoint^3)"]
    effectmatrix[x, 6] <- eff.fixed.timepoint3
    
    eff.random.intercept <- randomeffects[x, "(Intercept)"]
    effectmatrix[x, 7] <- eff.random.intercept 
    
    eff.random.timepoint <- randomeffects[x, "timepoint"]
    effectmatrix[x, 8] <- eff.random.timepoint 
    
    ind.marker <- genotypes[topmarker, x]
    eff.marker <- fixedeffects["mgt.AD"] * (as.numeric(ind.marker))
    if (is.na(ind.marker)) eff.marker <- NA
    effectmatrix[x, 9] <- eff.marker
    
    eff.markertime <- fixedeffects["timepoint:mgt.AD"] * (as.numeric(ind.marker))
    if (is.na(eff.markertime)) eff.markertime <- 0
    effectmatrix[x, 10] <- eff.markertime

    effectmatrix[x, 11] <- 0
    effectmatrix[x, 12] <- 0
    if(!is.na(ind.marker) && ind.marker != "0"){
      eff.markerDD <- fixedeffects["mgt.DD"]
      #if (is.na(ind.marker)) eff.markerDD <- NA
      effectmatrix[x, 11] <- eff.markerDD
      
      eff.markertimeDD <- fixedeffects["timepoint:mgt.DD"]
      if (is.na(eff.markertimeDD)) eff.markertimeDD <- 0
      effectmatrix[x, 12] <- eff.markertimeDD
    }
}
plot(c(0,120), c(0, 60), t = 'n')
estimates.null <- apply(effectmatrix, 1, function(eff){
  curve(eff["fIntercept"] + eff["Mother"] + eff["Sex"] + eff["tp"]*x + eff["tp2"]*x^2 + eff["tp3"]*x^3 + eff["marker"] + eff["markertime"]*x + eff["markerDD"] + eff["markertimeDD"]*x, from=0, to=120, add=TRUE)
})

plot(c(21,140), c(0, 60), t = 'n')
for(x in 1:nrow(phenotypes)){
  points(seq(21, 140, 7), phenotypes[x, paste0("d", seq(21, 140, 7))], t = 'b', col=rgb(0.1,0.1,0.1,0.1))
}
i <- 1
estimates.null <- apply(effectmatrix, 1, function(eff){
  a <- curve(eff["fIntercept"] + eff["Mother"] + eff["Sex"] + eff["tp"]*x + 
  eff["tp2"]*x^2 + eff["tp3"]*x^3 + eff["marker"] + eff["markertime"] * x + 
  eff["rIntercept"] + eff["markerDD"] + eff["markertimeDD"] * x, from=0, to=140, add=TRUE, col=as.numeric(genotypes[topmarker, rownames(effectmatrix)[i]]) + 2)
  i <<- i + 1
  return(a)
})



