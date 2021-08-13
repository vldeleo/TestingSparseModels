# Goal: Test sparse models against each other

  # things I would like to add but haven't yet:
    # BSLMM
    # RJ vs non-RJ
    # working calculations for monomvn's non-bayesian regression methods
    
### load modules ###
#module load swset/2018.05  gcc/7.3.0 r/4.0.5-py27
#py-matplotlib/2.2.2-py36 
#try using -m 50M with salloc?


### notes:
# regress() can do regression with partial least squares, principal component regression, lasso from the lars package, lar for least angle regression, ridge from the MASS package, and two classical forward covariate selection methods. cross validation in plsr and lars methods depend on random seed, but can use LOO instead.

# documentation for monomvn regress() SAYS it accepts dataframe or vector inputs, but it does not.

#######################################################
### load libraries ###
library(susieR)
library(monomvn)
library(glmnet)
library(dplyr)
library(reshape2)


### define modeling function ###
fitMod <- function(X, y, M=NULL, yOOS=NULL, XOOS=NULL, model = c("lm", "ols", "pls", "pcr", "lar", "lasso", "ridge", "blasso", "bridge", "bhs", "susie"), withPred=TRUE, coeff=FALSE){
  X <- X
  y <- y
  if(!is.null(M)){ M <- M } # can set max expected non-zero predictors; if sparsity is expected, lower values are faster. 

  
  # I saved my out of sample observations separately, so I will read that in and use
    #but I can create an test sample even without the pre-simulated set
if (withPred == TRUE){
    if (!is.null(yOOS) & !is.null(XOOS)) { # need both
    Xtest <- as.matrix(X)
    Xval <- as.matrix(XOOS)
    ytest <- as.matrix(y)
    yval <- as.matrix(yOOS)
  } else {   # create an out of sample 
    if (length(y)%%10 == 0){  
      valN <- length(y)/10   # I've split the sample 9:1, which works for any dataset with nObservations divisible by 10
    } else {valN <-  round(length(y)/10, 0)}   # for samples with a different length, I'll have to round 
    whichVal <- sample(1:length(y), valN)
    Xtest <- as.matrix(X[-(whichVal + 1),]) 
    Xval <- as.matrix(X[(whichVal + 1),])
    ytest <- as.matrix(y[-(whichVal)])
    yval <- as.matrix(y[whichVal])
    
    # check that the test and validation are similar to each other, resample if not
    # not sure if I want to include this step as a failsafe. need to think about it.
    while (t.test(ytest, yval)$p.value < 0.1) { # being conservative
      if (length(y)%%10 == 0){  
        valN <- length(y)/10   
      } else {valN <-  round(length(y)/10, 0)}    
      whichVal <- sample(1:length(y), valN)
      Xtest <- X[-(whichVal + 1),] 
      Xval <- X[(whichVal + 1),]
      ytest <- y[-(whichVal)]
      yval <- y[whichVal]
      }
} 
}
  
### fitting the models ###  
if ("lm" %in% model) {
  pt1 <- proc.time() # I'd like to compare how long each model takes
  if(exists("Xtest")){
    lmres <- lm(ytest ~ Xtest) } else {
      lmres <- lm(as.matrix(y) ~ as.matrix(X)) }
  lmTime <- as.numeric(proc.time() - pt1)[3]
  lmRMSE <- sqrt(mean((resid(lmres))^2))
  lmMAE <- mean(abs(resid(lmres)))
  lmAIC <- AIC(lmres)
  lmStats <- c(lmTime, lmRMSE, lmMAE, lmAIC)
}
if ("ols" %in% model) {
  pt1 <- proc.time() 
  if(exists("Xtest")){ # p = 0 forces regress to use my chosen method even when there are more columns than rows. can choose a value 0-1 to specify when to use lsr
    olsres <- regress(Xtest, ytest, method = c("lsr"), p = 0) 
    } else {
      olsres <- regress(as.matrix(X), as.matrix(y), method = c("lsr"), p = 0) }
  olsTime <- as.numeric(proc.time() - pt1)[3]
  olsRMSE <- sqrt(mean((resid(olsres))^2))
  olsMAE <- mean(abs(summary(olsres)$residuals))
  olsAIC <- AIC(olsres)
  olsStats <- c(olsTime, olsRMSE, olsMAE, olsAIC)
}
if ("pls" %in% model) {
  pt1 <- proc.time() 
  if(exists("Xtest")){ # can set max number of pc to consider
    plsres <- regress(Xtest, ytest, method = c("plsr"), p = 0, ncomp.max = 10)
    } else {
      plsres <- regress(as.matrix(X), as.matrix(y), method = c("plsr"), p = 0, ncomp.max = 10) }
  plsTime <- as.numeric(proc.time() - pt1)[3]
  plsRMSE <- sqrt(mean((summary(plsres)$residuals)^2))
  plsMAE <- mean(abs(summary(plsres)$residuals))
  plsAIC <- AIC(plsres)
  plsStats <- c(plsTime, plsRMSE, plsMAE, plsAIC)
}
if ("pcr" %in% model) {
  pt1 <- proc.time() 
  if(exists("Xtest")){ 
    pcrres <- regress(Xtest, ytest, method = c("pcr"), p = 0, ncomp.max = 10)
  } else {
    pcrres <- regress(as.matrix(X), as.matrix(y), method = c("pcr"), p = 0, ncomp.max = 10) }
    pcrTime <- as.numeric(proc.time() - pt1)[3]
    pcrRMSE <- sqrt(mean((summary(pcrres)$residuals)^2))
    pcrMAE <- mean(abs(summary(pcrres)$residuals))
    pcrAIC <- AIC(pcrres)
    pcrStats <- c(pcrTime, pcrRMSE, pcrMAE, pcrAIC)
}
if ("lar" %in% model) {
  pt1 <- proc.time() 
    if(exists("Xtest")){ 
    larres <- regress(Xtest, ytest, method = c("lar"), p = 0)
  } else {
    larres <- regress(as.matrix(X), as.matrix(y), method = c("lar"), p = 0) }
    larTime <- as.numeric(proc.time() - pt1)[3]
    larRMSE <- sqrt(mean((summary(larres)$residuals)^2))
    larMAE <- mean(abs(summary(larres)$residuals))
    larAIC <- AIC(larres)
    larStats <- c(larTime, larRMSE, larMAE, larAIC)
}
if ("lasso" %in% model) {
    pt1 <- proc.time() 
    if(exists("Xtest")){ 
      lasres <- regress(Xtest, ytest, method = c("lasso"), p = 0)
    } else {
      lasres <- regress(as.matrix(X), as.matrix(y), method = c("lasso"), p = 0) }
    lasTime <- as.numeric(proc.time() - pt1)[3]
    if(exists("Xtest")){lasResid <- c(lasres$mu + as.matrix(Xtest) %*% lasres$Beta) - ytest} else {
      lasResid <- mean(lasres$mu + as.matrix(X) %*% lasres$Beta) - y}
    lasRMSE <- sqrt(mean((lasResid)^2))
    lasMAE <- mean(abs(lasResid))
    lasAIC <- NA #AIC(lasres) I'm going to have to find a method for this
    lasStats <- c(lasTime, lasRMSE, lasMAE, lasAIC)
  }
if ("ridge" %in% model) {
    pt1 <- proc.time() 
    if(exists("Xtest")){ 
      ridgeres <- regress(Xtest, ytest, method = c("ridge"), p = 0)
    } else {ridgeres <- regress(as.matrix(X), as.matrix(y), method = c("ridge"), p = 0) }
    ridgeTime <- as.numeric(proc.time() - pt1)[3]
    ridgeRMSE <- sqrt(mean((summary(ridgeres)$residuals)^2))
    ridgeMAE <- mean(abs(summary(ridgeres)$residuals))
    ridgeAIC <- AIC(ridgeres)
    ridgeStats <- c(ridgeTime, ridgeRMSE, ridgeMAE, ridgeAIC)
  }
# getting bayesian #
if ("blasso" %in% model) {
  pt1 <- proc.time()
  if(exists("Xtest")){
    blassores <- blasso(Xtest, ytest, RJ=TRUE, M=M, T=25000, verb = 0)
  } else {blassores <- blasso(as.matrix(X), as.matrix(y), RJ=TRUE, M=M, T=25000, verb = 0)}
  blasTime <- as.numeric(proc.time() - pt1)[3]
  if(exists("Xtest")){blasResid <- c(mean(blassores$mu[10001:25000]) + as.matrix(Xtest) %*% colMeans(blassores$beta[10001:25000,])) - ytest} else {
    blasResid <- mean(blassores$mu[10001:25000]) + as.matrix(X) %*% colMeans(blassores$beta[10001:25000,]) - y}
  #blasResid <- pred(blassores(Xtest, ytest)) - ytest
  blasRMSE <- sqrt(mean((blasResid)^2))
  blasMAE <- mean(abs(blasResid))
  blasAIC <- NA #AIC(blassores) monomvn doesn't have a default function to do AIC on a blasso object
  blasStats <- c(blasTime, blasRMSE, blasMAE, blasAIC)
}
if ("bhs" %in% model) {
  pt1 <- proc.time()
  if(exists("Xtest")){bhsres <- bhs(Xtest, ytest, RJ=TRUE, M=M, T=25000, verb = 0)} else { 
    bhsres <- bhs(as.matrix(X), as.matrix(y), RJ=TRUE, M=M, T=25000, verb = 0)}
  bhsTime <- as.numeric(proc.time() - pt1)[3]
  if(exists("Xtest")){bhsResid <- c(mean(bhsres$mu[10001:25000]) + Xtest %*% colMeans(bhsres$beta[10001:25000,])) - ytest} else {
    bhsResid <- c(mean(bhsres$mu[10001:25000]) + as.matrix(X) %*% colMeans(bhsres$beta[10001:25000,])) - y}
  bhsRMSE <- sqrt(mean((bhsResid)^2))
  bhsMAE <- mean(abs(bhsResid))
  bhsAIC <- NA #AIC(bhsres) #need to find another approach here. is BIC comparable?
  bhsStats <- c(bhsTime, bhsRMSE, bhsMAE, bhsAIC)
}
if ("bridge" %in% model) {
  pt1 <- proc.time()
  if(exists("Xtest")){bridgeres <- bridge(Xtest, ytest, RJ=TRUE, M=M, T=25000, verb = 0)} else {bridgeres <- bridge(as.matrix(X), as.matrix(y), RJ=TRUE, M=M, T=25000, verb = 0)}
  briTime <- as.numeric(proc.time() - pt1)[3]
  if(exists("Xtest")){briResid <- c(mean(bridgeres$mu[10001:25000]) + Xtest %*% colMeans(bridgeres$beta[10001:25000,])) - ytest} else {
    briResid <- c(mean(bridgeres$mu[10001:25000]) + as.matrix(X) %*% colMeans(bridgeres$beta[10001:25000,])) - y}
  briRMSE <- sqrt(mean((briResid)^2))
  briMAE <- mean(abs(briResid))
  briAIC <- NA #AIC(brires)
  briStats <- c(briTime, briRMSE, briMAE, briAIC)
}
if ("susie" %in% model ) {
  pt1 <- proc.time()
  if(exists("Xtest")){
    susres <- susie(Xtest, ytest, L=M)} else { 
      susres <- susie(as.matrix(X), as.matrix(y), L=M) }
  susTime <- as.numeric(proc.time() - pt1)[3]
    if(exists("Xtest")){susResid <- predict(susres, type = "response") - ytest} else {susResid <- predict(susres, type = "response") - y}
  susRMSE <- sqrt(mean((susResid)^2))
  susMAE <- mean(abs(susResid))
  susAIC <- NA #AIC(susres) # looks like this one doesn't work either
  susStats <- c(susTime, susRMSE, susMAE, susAIC)
}
  

#predict y from models
if (withPred == TRUE) {
  if(exists("lmres")){
    lmpred <- lmres$coefficients[1] + Xval %*% lmres$coefficients[-1] 
    #note: R appears to still not have a good way to handle predictions of different shape than initial data, so I tried to do it the old fashioned way
    lmpredRMSE <- sqrt(mean((resid(lm(lmpred~yval)))^2))
    lmpredR2 <- summary(lm(lmpred~yval))$adj.r.squared
    lmpredcor <- cor.test(lmpred, yval)
    lmPStats <- c(lmpredRMSE, lmpredR2, lmpredcor)
  }
  if(exists("olsres")){
    olspred <- olsres$b[1] + Xval %*% olsres$b[-1]
    olspredRMSE <- sqrt(mean((summary(lm(olspred~yval))$residuals)^2))
    olspredR2 <- summary(lm(olspred~yval))$adj.r.squared
    olspredcor <- cor.test(olspred, yval) # are lm() and regress(method = "ols") equivalent? it looks like it
    olsPStats <- c(olspredRMSE, olspredR2, olspredcor)
  }
  if(exists("plsres")){
    plspred <- plsres$b[1] + Xval %*% plsres$b[-1]
    plspredRMSE <- sqrt(mean((summary(lm(plspred~yval))$residuals)^2))
    plspredR2 <- summary(lm(plspred~yval))$adj.r.squared
    plspredcor <- cor.test(plspred, yval) 
    plsPStats <- c(plspredRMSE, plspredR2, plspredcor)
  }
  if(exists("pcrres")){
    pcrpred <- pcrres$b[1] + Xval %*% pcrres$b[-1]
    pcrpredRMSE <- sqrt(mean((summary(lm(pcrpred~yval))$residuals)^2))
    pcrpredR2 <- summary(lm(pcrpred~yval))$adj.r.squared
    pcrpredcor <- cor.test(pcrpred, yval)
    pcrPStats <- c(pcrpredRMSE, pcrpredR2, pcrpredcor)
  }
  if(exists("larres")){
    larpred <- larres$b[1] + Xval %*% larres$b[-1]
    larpredRMSE <- sqrt(mean((summary(lm(larpred~yval))$residuals)^2))
    larpredR2 <- summary(lm(larpred~yval))$adj.r.squared
    larpredcor <- cor.test(larpred, yval) 
    larPStats <- c(larpredRMSE, larpredR2, larpredcor)
  }
  if(exists("lasres")){
    lassopred <- coef(lassores)[1] + Xval %*% coef(lassores)[-1]
    lassopredRMSE <- sqrt(mean((summary(lm(lassopred~yval))$residuals)^2))
    lassopredR2 <- summary(lm(lassopred~yval))$adj.r.squared
    lassopredcor <- cor.test(lassopred, yval) 
    lassoPStats <- c(lassopredRMSE, lassopredR2, lassopredcor)
  }
  if(exists("ridgeres")){
    ridgepred <- coef(ridgeres)[1] + Xval %*% coef(ridgeres)[-1]
    ridgepredRMSE <- sqrt(mean((summary(lm(ridgepred~yval))$residuals)^2))
    ridgepredR2 <- summary(lm(ridgepred~yval))$adj.r.squared
    ridgepredcor <- cor.test(ridgepred, yval) 
    ridgePStats <- c(ridgepredRMSE, ridgepredR2, ridgepredcor)
  }
  if(exists("blassores")){
    blassopred <- mean(blassores$mu[10001:25000]) + Xval %*% colMeans(blassores$beta[10001:25000,])
    blaspredRMSE <- sqrt(mean((resid(lm(blassopred~yval)))^2))
    blaspredR2 <- summary(lm(blassopred~yval))$adj.r.squared
    blaspredcor <- cor.test(blassopred, yval) 
    blasPStats <- c(blaspredRMSE, blaspredR2, blaspredcor)
  }
  if(exists("bridgeres")){
    bridgepred <- mean(bridgeres$mu[10001:25000]) + Xval %*% colMeans(bridgeres$beta[10001:25000,]) 
    bripredRMSE <- sqrt(mean((summary(lm(bridgepred~yval))$residuals)^2))
    bripredR2 <- summary(lm(bridgepred~yval))$adj.r.squared
    bripredcor <- cor.test(bridgepred, yval) 
    briPStats <- c(bripredRMSE, bripredR2, bripredcor)
  }
  if(exists("bhsres")){
    bhspred <- mean(bhsres$mu[10001:25000]) + Xval %*% colMeans(bhsres$beta[10001:25000,]) 
    bhspredRMSE <- sqrt(mean((summary(lm(bhspred~yval))$residuals)^2))
    bhspredR2 <- summary(lm(bhspred~yval))$adj.r.squared
    bhspredcor <- cor.test(bhspred, yval) 
    bhsPStats <- c(bhspredRMSE, bhspredR2, bhspredcor)
  }
  if(exists("susres")){
    suspred <- coef(susres)[1] + Xval %*% coef(susres)[-1]
    #going back to my note above, predict doesn't work with a differently shaped NewData. I tried doubling my oos so it would be the same size, but it looks like pred still used x values from my original set.
    suspredRMSE <- sqrt(mean((summary(lm(suspred~yval))$residuals)^2))
    suspredR2 <- summary(lm(suspred~yval))$adj.r.squared
    suspredcor <- cor.test(suspred, yval) 
    susPStats <- c(suspredRMSE, suspredR2, suspredcor)
  }
}
  
# it's quicker and more efficient to output performance metrics rather than model parameters, but I'll leave that as an option anyway
  
#get output
input_out <- list(if(exists("ytest")){ytest}, 
                  if(exists("yval")){yval}, 
                  if(exists("Xtest")){Xtest}, 
                  if(exists("Xval")){Xval})
names(input_out) <- c("ytest", "yval", "Xtest", "Xval")
input_out <- input_out[!sapply(input_out,is.null)]

stats_out <- list(if(exists("lmres")){lmStats}, 
                  if(exists("olsres")){olsStats}, 
                  if(exists("plsres")){plsStats}, 
                  if(exists("pcrres")){pcrStats}, 
                  if(exists("larres")){larStats}, 
                  if(exists("lasres")){lasStats}, 
                  if(exists("ridgeres")){ridgeStats}, 
                  if(exists("blassores")){blasStats}, 
                  if(exists("bridgeres")){briStats}, 
                  if(exists("bhsres")){bhsStats}, 
                  if(exists("susres")){susStats})
names(stats_out) <- c("lmStats", "olsStats", "plsStats", "pcrStats", "larStats", "lassoStats", "ridgeStats", "blassoStats", "bridgeStats", "bhsStats", "susStats")
stats_out <- stats_out[!sapply(stats_out,is.null)]

if (withPred==TRUE){
  pred_out <- list(if(exists("lmPStats")){lmPStats}, 
                   if(exists("olsPStats")){olsPStats}, 
                   if(exists("plsPStats")){plsPStats}, 
                   if(exists("pcrPStats")){pcrPStats}, 
                   if(exists("larPStats")){larPStats}, 
                   if(exists("lassoPStats")){lassoPStats}, 
                   if(exists("ridgePStats")){ridgePStats}, 
                   if(exists("blassoPStats")){blassoPStats}, 
                   if(exists("bridgePStats")){bridgePStats}, 
                   if(exists("bhsPStats")){bhsPStats}, 
                   if(exists("susPStats")){susPStats})
  names(pred_out) <- c("lmPStats", "olsPStats", "plsPStats", "pcrPStats", "larPStats", "lassoPStats", "ridgePStats", "blasPStats", "bridgePStats", "bhsPStats", "susPStats")
  pred_out <- pred_out[!sapply(pred_out,is.null)]
}

if (coeff==TRUE) { # IF user wants to keep all the estimated coefficients, this is how they could do that
coeff_out <- list(if(exists("lmres")){coef(lmres)}, 
                  if(exists("olsres")){olsres$b}, 
                  if(exists("plsres")){plsres$b}, 
                  if(exists("pcrres")){pcrres$b}, 
                  if(exists("larres")){larres$b}, 
                  if(exists("lassores")){coef(lassores)}, 
                  if(exists("ridgeres")){coef(ridgeres)}, 
                  if(exists("blassores")){c(mean(blassores$mu[10001:25000]),colMeans(blassores$beta[10001:25000,]))}, 
                  if(exists("bridgeres")){c(mean(bridgeres$mu[10001:25000]),colMeans(bridgeres$beta[10001:25000,]))}, 
                  if(exists("bhsres")){c(mean(bhsres$mu[10001:25000]),colMeans(bhsres$beta[10001:25000,]))}, 
                  if(exists("susres")){coef(susres)})
names(coeff_out) <- c("lmBeta", "olsBeta", "plsBeta", "pcrBeta", "larBeta", "lassoBeta", "ridgeBeta", "blasBeta", "bridgeBeta", "bhsBeta", "susBeta")
coeff_out <- coeff_out[!sapply(coeff_out,is.null)]
}

#put it together (but do I want it together?)
out_alltog <- list(input_out,
                   stats_out,
                   if(exists("coeff_out")){coeff_out},
                   if(exists("pred_out")){pred_out})
names(out_alltog) <- c("Input", "ModelStats", "Coefficients", "PredictionStats")
out_alltog <- out_alltog[!sapply(out_alltog, is.null)]

return(out_alltog)  # potentially a huge, confusing output
}




###################################################################
### now run it ###
### rewrote this so that I don't need to load all my simulated data into memory at the same time and crash my session ###
resFiles <- c(list.files(pattern = "simy.csv", full.names=TRUE)) # create list of files containing my response variables
predFiles <- c(list.files(pattern = "simX.csv", full.names=TRUE)) # create list of files containing my predictor variables

# check that there are an equal number of predictors and responses files
if (length(resFiles) != length(predFiles)){
  print("Predictors and variables are not the same length. Check input.")
}

#*** Optional, if I want to use an already sample of extra data ***
yOOSFiles <- c(list.files(pattern = "yOOS.csv", full.names=TRUE)) # create a list of files containing my reserved out of sample y values
XOOSFiles <- c(list.files(pattern = "XOOS.csv", full.names=TRUE)) 


# fit models over the replicates of simulated data in my directory
for (l in 1:length(resFiles)){
  # load data, 1 rep at a time
  resData <- as.data.frame(read.csv(file=resFiles[l], header = T))[,-1, drop = FALSE]
  predData <- read.csv(file=predFiles[l],header=T)[,-1]
  # check that predictor datasets match length of response datasets; if not, send warning
  if (nrow(resData) != nrow(predData)){
    print("Predictors and variables do not match. Check input.")
  }
  
  # if I want to use my already simulated out of sample observations
  if(exists("yOOSFiles")){
    yOOS <- as.data.frame(read.csv(file=yOOSFiles[l],header=T))[,-1, drop = FALSE]
  }
  if(exists("XOOSFiles")){
    XOOS <- read.csv(file=XOOSFiles[l],header=T)[,-1]
  }
       
  # fit the models, make predictions, return model diagnostics and coefficients (if requested)
  modOut <- fitMod(predData, resData, 15, yOOS, XOOS, model = c("lm", "pls", "pcr", "lar", "lasso", "ridge", "blasso", "bridge", "bhs", "susie"), withPred=TRUE, coeff=FALSE) 
  
  # save output
  for (c in 1:length(modOut)){ 
    write.csv(modOut[c], file = paste("ModelResults/Rep", l, "_", names(modOut)[c], ".csv", sep = ""))
  }
  #rm(list=ls()) looking for a way to clean up intermediary files so my process doesn't get killed... this might be too much, though, since it deletes my function
  gc() #this is not enough. Still get a segmentation fault with > 1 rep
  print(l)
}

#################################################



############ Scratch ##########################
# glmnet, and generalized sparse linear models
if ("glasso" %in% model){
  pt1 <- proc.time()
  if(exists(xtest)){lassores <- glmnet(Xtest[,-1], ytest, alpha=1,verb = 0)} else {lassores <- glmnet(X[,-1], y, verb = 0)}
  lasTime <- proc.time() - pt1
}
if ("gridge" %in% model){
  pt1 <- proc.time()
  if(exists(xtest)){ridgeres <- glmnet(Xtest[,-1], ytest, alpha=0, verb = 0)} else {ridgeres <- glmnet(X[,-1], y, RJ=FALSE, M=M, T=25000, verb = 0)}
  ridgeTime <- proc.time() - pt1
}
if ("elnet" %in% model){
  pt1 <- proc.time()
  if(exists(xtest)){elnetres <- glmnet(Xtest[,-1], ytest, alpha=.5, verb = 0)} else {elnetres <- glmnet(X[,-1], y, verb = 0)} #I can specify the elasticnet mixing parameter between 1 and 0; I've defaulted to midway between
  elnetTime <- proc.time() - pt1
}
#I don't fully understand glmnet's "relax" option, so I've left it out. documentation says that this feature is best if fitting the original glmnet model takes a long time




