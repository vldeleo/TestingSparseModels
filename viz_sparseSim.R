# Now, I can look at some statistics and compare performance of the models I made

# Because I have to remind myself of this line EVERYTIME I log on to the server:
  #module load swset/2018.05  gcc/7.3.0 py-matplotlib/2.2.2-py36 r/4.0.5-py27

# for each rep 1-100, I have "simy", "simX", "simBeta", "yOOS", "XOOS", "lmBeta", "lmPred", "susBeta", "susPred", "blasBeta", "blasPred", "bhsBeta", "bhsPred", "bridgeBeta", "bridgePred"

# set up R
library(rethinking)
library(dplyr)
library(reshape2)

# compile my data for Betas
simBeta_fn <- list.files(pattern = "simBeta.csv", full.names=TRUE)
Betafew_data = lapply(simBeta_fn, function(x){read.csv(file=x,header=T)})
simBetaall <- Reduce(function(x,y){full_join(x,y, by = "X")}, Betafew_data)
colnames(simBetaall) <- c("X", paste("rep",1:100, sep = ""))
#
lmBeta_fn <- list.files(pattern = "lmBeta.csv", full.names=TRUE)
lmBeta_data = lapply(lmBeta_fn, function(x){read.csv(file=x,header=T)})
lmBetaall <- Reduce(function(x,y){full_join(x,y, by = "X")}, lmBeta_data)
colnames(lmBetaall) <- c("X", paste("rep",1:100, sep = ""))
#
susBeta_fn <- list.files(pattern = "susBeta.csv", full.names=TRUE)
susBeta_data = lapply(susBeta_fn, function(x){read.csv(file=x,header=T)})
susBetaall <- Reduce(function(x,y){full_join(x,y, by = "X")}, susBeta_data)
colnames(susBetaall) <- c("X", paste("rep",1:100, sep = ""))
#
blasBeta_fn <- list.files(pattern = "blasBeta.csv", full.names=TRUE)
blasBeta_data = lapply(blasBeta_fn, function(x){read.csv(file=x,header=T)})
blasBetaall <- Reduce(function(x,y){full_join(x,y, by = "X")}, blasBeta_data)
colnames(blasBetaall) <- c("X", paste("rep",1:100, sep = ""))
#
bhsBeta_fn <- list.files(pattern = "bhsBeta.csv", full.names=TRUE)
bhsBeta_data = lapply(bhsBeta_fn, function(x){read.csv(file=x,header=T)})
bhsBetaall <- Reduce(function(x,y){full_join(x,y, by = "X")}, bhsBeta_data)
colnames(bhsBetaall) <- c("X", paste("rep",1:100, sep = ""))
#
briBeta_fn <- list.files(pattern = "bridgeBeta.csv", full.names=TRUE)
briBeta_data = lapply(briBeta_fn, function(x){read.csv(file=x,header=T)})
briBetaall <- Reduce(function(x,y){full_join(x,y, by = "X")}, briBeta_data)
colnames(briBetaall) <- c("X", paste("rep",1:100, sep = ""))
#


### pull up the data used for OOS predictions ###
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
#
susPredfew_fn <- list.files(pattern = "susPred.csv", full.names=TRUE) 
susPredfew_data = lapply(susPredfew_fn, function(x){read.csv(file=x,header=T)})
susPredall <- Reduce(function(x,y){full_join(x,y, by = "X")}, susPredfew_data)
colnames(susPredall) <- c("X", paste("rep",1:100, sep = ""))
#
blasPredfew_fn <- list.files(pattern = "blasPred.csv", full.names=TRUE) 
blasPredfew_data = lapply(blasPredfew_fn, function(x){read.csv(file=x,header=T)})
blasPredall <- Reduce(function(x,y){full_join(x,y, by = "X")}, blasPredfew_data)
colnames(blasPredall) <- c("X", paste("rep",1:100, sep = ""))
#
bhsPredfew_fn <- list.files(pattern = "bhsPred.csv", full.names=TRUE) 
bhsPredfew_data = lapply(bhsPredfew_fn, function(x){read.csv(file=x,header=T)})
bhsPredall <- Reduce(function(x,y){full_join(x,y, by = "X")}, bhsPredfew_data)
colnames(bhsPredall) <- c("X", paste("rep",1:100, sep = ""))
#
briPredfew_fn <- list.files(pattern = "bridgePred.csv", full.names=TRUE) 
briPredfew_data = lapply(briPredfew_fn, function(x){read.csv(file=x,header=T)})
briPredall <- Reduce(function(x,y){full_join(x,y, by = "X")}, briPredfew_data)
colnames(briPredall) <- c("X", paste("rep",1:100, sep = ""))





# 1. How well do model coefficients match true coefficients?
# Pearsons?
  # I'm not sure my data are linear + homoscedastic (actually, almost positive these won't be homoscedastic for the dense simulated model)
  # if I ranked my effect sizes, I could try Spearman's coefficient- I do hope for a monotonic relationship here.
shapiroReal <- lapply(2:101,function(rep){shapiro.test(simBetaall[,rep])})
shapiroLM <- lapply(2:101,function(rep){shapiro.test(lmBetaall[2:401,rep])})
      # hmmm very non-normal

lmBetaRho <- unlist(lapply(2:101,function(rep){unname(cor.test(simBetaall[,rep], lmBetaall[2:401,rep])$estimate)}))
susBetaRho <- unlist(lapply(2:101,function(rep){unname(cor.test(simBetaall[,rep], susBetaall[2:401,rep])$estimate)}))
blasBetaRho <- unlist(lapply(2:101,function(rep){unname(cor.test(simBetaall[,rep], blasBetaall[2:401,rep])$estimate)}))
bhsBetaRho <- unlist(lapply(2:101,function(rep){unname(cor.test(simBetaall[,rep], bhsBetaall[2:401,rep])$estimate)}))
briBetaRho <- unlist(lapply(2:101,function(rep){unname(cor.test(simBetaall[,rep], briBetaall[2:401,rep])$estimate)}))

BetaRhotog <- as.data.frame(cbind(lmBetaRho, susBetaRho, blasBetaRho, bhsBetaRho, briBetaRho))
BetaRholong <- melt(BetaRhotog, variable.name = "Model")

png(filename = "many_BetaRhoboxplot_5.28.21.png", width = 600, height = 400)
boxplot(value ~ Model, data = BetaRholong, main = "")
dev.off()

# I wrote the model vs actual df to file so I could upload to github for Rmarkdown
write.csv(BetaRhotog, "many_BetaRho.csv")


# Root Mean Squared Error of Betas
lmBetaRMSE <- unlist(lapply( 2:101, function(rep){sqrt(mean((summary(lm(lmBetaall[2:401,rep]~simBetaall[,rep]))$residuals)^2))}))
susBetaRMSE <- unlist(lapply( 2:101, function(rep){sqrt(mean((summary(lm(susBetaall[2:401,rep]~simBetaall[,rep]))$residuals)^2))}))
blasBetaRMSE <- unlist(lapply( 2:101, function(rep){sqrt(mean((summary(lm(blasBetaall[2:401,rep]~simBetaall[,rep]))$residuals)^2))}))
bhsBetaRMSE <- unlist(lapply( 2:101, function(rep){sqrt(mean((summary(lm(bhsBetaall[2:401,rep]~simBetaall[,rep]))$residuals)^2))}))
briBetaRMSE <- unlist(lapply( 2:101, function(rep){sqrt(mean((summary(lm(briBetaall[2:401,rep]~simBetaall[,rep]))$residuals)^2))}))
# get a warning about the bridge predictions being "essentially perfect fit: summary may be unreliable"
BetaRMSEtog <- as.data.frame(cbind(lmBetaRMSE, susBetaRMSE, blasBetaRMSE, bhsBetaRMSE, briBetaRMSE))
BetaRMSElong <- melt(BetaRMSEtog, variable.name = "Model")

png(filename = "many_BetaRMSEboxplot_5.28.21.png", width = 600, height = 400)
boxplot(value ~ Model, data = BetaRMSElong)
dev.off()

# writing the file
write.csv(BetaRMSEtog, "many_BetaRMSE.csv")


#r^2
lmBetaR2 <- unlist(lapply( 2:101, function(rep){summary(lm(lmBetaall[2:401,rep]~simBetaall[,rep]))$adj.r.squared})) # be sure to skip first row of model predictions, which should be the intercept
susBetaR2 <- unlist(lapply( 1:100, function(rep){summary(lm(susBetaall[2:401,rep+1]~simBetaall[,rep+1]))$adj.r.squared}))
blasBetaR2 <- unlist(lapply( 1:100, function(rep){summary(lm(blasBetaall[2:401,rep+1]~simBetaall[,rep+1]))$adj.r.squared}))
bhsBetaR2 <- unlist(lapply( 1:100, function(rep){summary(lm(bhsBetaall[2:401,rep+1]~simBetaall[,rep+1]))$adj.r.squared}))
briBetaR2 <- unlist(lapply( 1:100, function(rep){summary(lm(briBetaall[2:401,rep+1]~simBetaall[,rep+1]))$adj.r.squared}))
# 
BetaR2tog <- as.data.frame(cbind(lmBetaR2, susBetaR2, blasBetaR2, bhsBetaR2, briBetaR2))
BetaR2long <- melt(BetaR2tog, variable.name = "Model")

png(filename = "many_BetaR2boxplot_5.28.21.png", width = 600, height = 400)
boxplot(value ~ Model, data = BetaR2long, main = "")
dev.off()

# I wrote the model vs actual df to file so I could upload to github for Rmarkdown
write.csv(BetaR2tog, "many_BetaR2.csv")

  # for the models fit to data with many true effects, the r2 may be misleading because off course many effects are set to 0. would be interesting to see how this differs at different effect sizes.
lmBetaHE <- lmBetaall[lmBetaall > (mean(lmBetall) + sd(lmBetaall)) | lmBetaall < (mean(lmBetall) + sd(lmBetaall))]

lmBetaR2 <- unlist(lapply( 1:100, function(rep){summary(lm(lmBetaall[2:401,rep+1]~simBetaall[,rep+1]))$adj.r.squared})) # be sure to skip first row of model predictions, which should be the intercept
susBetaR2 <- unlist(lapply( 1:100, function(rep){summary(lm(susBetaall[2:401,rep+1]~simBetaall[,rep+1]))$adj.r.squared}))
blasBetaR2 <- unlist(lapply( 1:100, function(rep){summary(lm(blasBetaall[2:401,rep+1]~simBetaall[,rep+1]))$adj.r.squared}))
bhsBetaR2 <- unlist(lapply( 1:100, function(rep){summary(lm(bhsBetaall[2:401,rep+1]~simBetaall[,rep+1]))$adj.r.squared}))
briBetaR2 <- unlist(lapply( 1:100, function(rep){summary(lm(briBetaall[2:401,rep+1]~simBetaall[,rep+1]))$adj.r.squared}))
# 
BetaR2tog <- as.data.frame(cbind(lmBetaR2, susBetaR2, blasBetaR2, bhsBetaR2, briBetaR2))
BetaR2long <- melt(BetaR2tog, variable.name = "Model")




# 2. How well do model predictions match OOS true response variables?


# Root Mean Squared Error of predictions
lmPredRMSE <- unlist(lapply( 2:101, function(rep){sqrt(mean((summary(lm(lmPredall[,rep]~yOOSall[,rep]))$residuals)^2))}))
susPredRMSE <- unlist(lapply( 2:101, function(rep){sqrt(mean((summary(lm(susPredall[,rep]~yOOSall[,rep]))$residuals)^2))}))
susPredRMSE <- unlist(lapply( 2:101, function(rep){sqrt(mean((summary(lm(susPredall[,rep]~yOOSall[,rep]))$residuals)^2))}))
blasPredRMSE <- unlist(lapply( 2:101, function(rep){sqrt(mean((summary(lm(blasPredall[,rep]~yOOSall[,rep]))$residuals)^2))}))
bhsPredRMSE <- unlist(lapply( 2:101, function(rep){sqrt(mean((summary(lm(bhsPredall[,rep]~yOOSall[,rep]))$residuals)^2))}))
briPredRMSE <- unlist(lapply( 2:101, function(rep){sqrt(mean((summary(lm(briPredall[,rep]~yOOSall[,rep]))$residuals)^2))}))
# get a warning about the bridge predictions being "essentially perfect fit: summary may be unreliable"
predRMSEtog <- as.data.frame(cbind(lmPredRMSE, susPredRMSE, blasPredRMSE, bhsPredRMSE, briPredRMSE))
predRMSElong <- melt(predRMSEtog, variable.name = "Model")

png(filename = "many_PredRMSEboxplot_5.28.21.png", width = 600, height = 400)
boxplot(value ~ Model, data = predRMSElong)
dev.off()

# writing the file
write.csv(predRMSEtog, "many_PredRMSE.csv")



# using R2
lmPredR2 <- unlist(lapply( 1:100, function(rep){summary(lm(lmPredall[,rep+1]~yOOSall[,rep+1]))$adj.r.squared}))
susPredR2 <- unlist(lapply( 1:100, function(rep){summary(lm(susPredall[,rep+1]~yOOSall[,rep+1]))$adj.r.squared}))
susPredR2 <- unlist(lapply( 1:100, function(rep){summary(lm(susPredall[,rep+1]~yOOSall[,rep+1]))$adj.r.squared}))
blasPredR2 <- unlist(lapply( 1:100, function(rep){summary(lm(blasPredall[,rep+1]~yOOSall[,rep+1]))$adj.r.squared}))
bhsPredR2 <- unlist(lapply( 1:100, function(rep){summary(lm(bhsPredall[,rep+1]~yOOSall[,rep+1]))$adj.r.squared}))
briPredR2 <- unlist(lapply( 1:100, function(rep){summary(lm(briPredall[,rep+1]~yOOSall[,rep+1]))$adj.r.squared}))
  # get a warning about the bridge predictions being "essentially perfect fit: summary may be unreliable"
predR2tog <- as.data.frame(cbind(lmPredR2, susPredR2, blasPredR2, bhsPredR2, briPredR2))
predR2long <- melt(predR2tog, variable.name = "Model")

png(filename = "many_R2boxplot_5.28.21.png", width = 600, height = 400)
boxplot(value ~ Model, data = predR2long)
dev.off()

# I wrote the model vs prediction df to file so I could upload to github for Rmarkdown
write.csv(predR2tog, "many_PredR2.csv")



# 3. really basic, what are the performance stats on the models?
#Root Mean Squared Error:
d = simy-pred(lmBeta(simX))
mse = mean((d)^2)
mae = mean(abs(d))
rmse = sqrt(mse)


# in reverse:
rmse <- sqrt(mean((simy-)^2))

sqrt(mean(M$residuals^2))

