# Now, I can look at some statistics and compare performance of the models I made

# Because I have to remind myself of this line EVERYTIME I log on to the server:
  #module load swset/2018.05  gcc/7.3.0 py-matplotlib/2.2.2-py36 r/4.0.5-py27

# for each rep 1-100, I have "simy", "simX", "simBeta", "yOOS", "XOOS", "lmBeta", "lmPred", "susBeta", "susPred", "blasBeta", "blasPred", "bhsBeta", "bhsPred", "bridgeBeta", "bridgePred"

# set up R
library(rethinking)
library(dplyr)

# compile my data
yfew_files <- list.files(pattern = "simy.csv", full.names=TRUE) 
#test <- read.csv(allfiles[1], header = TRUE)
yfew_list = lapply(yfew_files, function(x){read.csv(file=x,header=T)})
Reduce(function(x,y) {merge(x,y, by = "row.names")}, yfew_list)
  # ^ instead of this, use the code I figured out below

#
yOOSfew_fn <- list.files(pattern = "yOOS.csv", full.names=TRUE) 
yOOSfew_data <- lapply(yOOSfew_fn, function(x){read.csv(file=x,header=T)})
yOOSall <- Reduce(function(x,y){full_join(x,y, by = "X")}, yOOSfew_data)
colnames(yOOSall) <- c("X", paste("rep",1:100, sep = ""))

#
lmPredfew_fn <- list.files(pattern = "lmPred.csv", full.names=TRUE) 
lmPredfew_data = lapply(lmPredfew_fn, function(x){read.csv(file=x,header=T)})
lmPredall <- Reduce(function(x,y){full_join(x,y, by = "X")}, lmPredfew_data)
colnames(lmPredall) <- c("X", paste("rep",1:100, sep = ""))


# 1. How well do model coefficients match true coefficients?

lmCoefCor <- apply(function(rep){cor.test(simBeta[rep], lmBeta[rep])}, 1:100)




# 2. How well do model predictions match OOS true response variables?

lmPredCor <- apply(function(rep){cor.test(yOOS[rep], lmPred[rep])}, 1:100)
susPredCor <- apply(function(rep){cor.test(yOOS[rep], susPred[rep])}, 1:100)


boxplot(lmPredCor, susPredCor)


# using R2 (preferred method)
lmPredR2 <- lapply( 1:100, function(rep){summary(lm(lmPredall[,rep+1]~yOOSall[,rep+1]))$adj.r.squared})


# Brier score? (for classification problems, not relevant)
