setwd("D:/Edrive/Mouse/AIL_Manuel")

phenotypes <- read.csv("phenotypes.txt", header=TRUE, check.names=FALSE, sep="\t", row.names = 1)
measurementdays.char <- paste0("d", seq(21, 140, 7))
measurementdays.num <- as.numeric(gsub("d", "", measurementdays.char))

weight <- phenotypes[, measurementdays.char]

#QC by weight
plot(x = c(min(measurementdays.num), max(measurementdays.num)), 
     y = c(0, max(weight * 1.15, na.rm=TRUE)), 
     main = "Individual growth curves",
     t = "n", yaxs = "i", xlab = "time(days)", ylab = "weight(gr.)")
for (x in 1:nrow(phenotypes)) {
  points(measurementdays.num, weight[x,], t ="l", col = "orange")
}

#QC by weight
plot(x = c(100, 140), 
     y = c(0, max(weight * 1.15, na.rm=TRUE)), 
     main = "Individual growth curves",
     t = "n", yaxs = "i", xlab = "time(days)", ylab = "weight(gr.)")
for (x in 1:nrow(phenotypes)) {
  points(measurementdays.num, weight[x,], t ="l", col = "orange")
}

# oralGTT
oralGTT.cols <- paste0("GTT_", c(0, 15, 30, 60, 120))
oralGTT <- phenotypes[, oralGTT.cols]
oralGTT.num <- as.numeric(gsub("GTT_", "", colnames(oralGTT)))

plot(main ="Oral Glucose Tolerance Test (GTT)", 
     x = c(min(oralGTT.num), 
     y = max(oralGTT.num)), c(0, max(oralGTT) * 1.15),  t = 'n', xlab="Time (min)", ylab="Blood Glucose (mg/dl)")

for (x in 1:nrow(oralGTT)) {
  points(oralGTT.num, oralGTT[x,], t="l", col="orange")
}

oralGTT.boxplot <- boxplot(oralGTT, main="Oral Glucose Tolerance Test", ylab="Blood Glucose (mg/dl)", xlab="Time(min)", col = "orange", notch=TRUE)
lines(1:5, oralGTT.boxplot$stats[ 3, ], col="blue", lwd=2)

# insulinTT
ITT.cols <-  paste0("ITT_", c(0,15,30,60))
ITT <- phenotypes[, ITT.cols]
ITT.num <- as.numeric(gsub("ITT_", "", colnames(ITT)))

plot(main ="Insulin Tolerance Test (ITT)", 
     x = c(min(ITT.num), max(ITT.num)), 
     y = c(0, max(ITT) * 1.15),  t = 'n', xlab="Time (min)", ylab="Blood Glucose (mg/dl)")

for (x in 1:nrow(ITT)){
  points(ITT.num, ITT[x,], t="l", col="orange")
}

ITT.boxplot <- boxplot(ITT, main="Insulin Tolerance Test", ylab="Blood Glucose (mg/dl)", xlab="Time(min)", col = "orange", notch=TRUE)
lines(1:4, ITT.boxplot$stats[ 3, ], col="blue", lwd=2)

#Tissues weight
tissues <- c("Leber", "WATgon", "WATsc")
tissues.weights <- phenotypes[, tissues]
ordering <- sort(tissues.weights[,"WATgon"], index.return=TRUE)$ix
tissues.weights <- tissues.weights[ordering,]

plot(main="weight relationship", c(0, nrow(tissues.weights)), c(0, max(tissues.weights, na.rm=TRUE)))
lines(tissues.weights[,"WATgon"], col = "red" , lwd=3 , pch=19 , type="l")
lines(tissues.weights[,"Leber"], col = "blue" , lwd=3 , pch=19 , type="l")
lines(tissues.weights[,"WATsc"], col = "green" , lwd=3 , pch=19 , type="l")

legend("topleft",
  legend = c("Liver", "Gon", "SCF"),
  col = c("blue", "red", "green"),
  pch = c(20,20,20),
  bty = "n",
  pt.cex = 2,
  cex = 1.2,
  text.col = "black")
  