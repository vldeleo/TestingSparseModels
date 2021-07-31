# Test sparse models against each other

  # things I would like to add but haven't yet:
    # BSLMM
    # RJ vs non-RJ

### load modules ###
#module load swset/2018.05  gcc/7.3.0 r/4.0.5-py27
#py-matplotlib/2.2.2-py36 
#try using -m 50M with salloc?


# regress() can do regression with partial least squares, principal component regression, lasso from the lars package, lar for least angle regression, ridge from the MASS package, and two classical forward covariate selection methods. cross validation in plsr and lars methods depend on random seed, but can use LOO instead.

# documentation for monomvn regress() SAYS it accepts dataframe or vector inputs, but it does not.


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
  lmRMSE <- 
  lmMAE <- 
  lmAIC <- 
  lmStats <- c(lmRMSE, lmMAE, lm)
}
if ("ols" %in% model) {
  pt1 <- proc.time() 
  if(exists("Xtest")){ # p = 0 forces regress to use my chosen method even when there are more columns than rows. can choose a value 0-1 to specify when to use lsr
    olsres <- regress(Xtest, ytest, method = c("lsr"), p = 0) 
    } else {
      olsres <- regress(as.matrix(X), as.matrix(y), method = c("lsr"), p = 0) }
  olsTime <- as.numeric(proc.time() - pt1)[3]
}
if ("pls" %in% model) {
  pt1 <- proc.time() 
  if(exists("Xtest")){ # can set max number of pc to consider
    plsres <- regress(Xtest, ytest, method = c("plsr"), p = 0, ncomp.max = 10)
    } else {
      plsres <- regress(as.matrix(X), as.matrix(y), method = c("plsr"), p = 0, ncomp.max = 10) }
  plsTime <- as.numeric(proc.time() - pt1)[3]
}
if ("pcr" %in% model) {
  pt1 <- proc.time() 
  if(exists("Xtest")){ 
    pcrres <- regress(Xtest, ytest, method = c("pcr"), p = 0, ncomp.max = 10)
  } else {
    pcrres <- regress(as.matrix(X), as.matrix(y), method = c("pcr"), p = 0, ncomp.max = 10) }
    pcrTime <- as.numeric(proc.time() - pt1)[3]
}
if ("lar" %in% model) {
  pt1 <- proc.time() 
    if(exists("Xtest")){ 
    larres <- regress(Xtest, ytest, method = c("lar"), p = 0)
  } else {
    larres <- regress(as.matrix(X), as.matrix(y), method = c("lar"), p = 0) }
    larTime <- as.numeric(proc.time() - pt1)[3]
}
if ("lasso" %in% model) {
    pt1 <- proc.time() 
    if(exists("Xtest")){ 
      lasres <- regress(Xtest, ytest, method = c("lasso"), p = 0)
    } else {
      lasres <- regress(as.matrix(X), as.matrix(y), method = c("lasso"), p = 0) }
    lasTime <- as.numeric(proc.time() - pt1)[3]
  }
if ("ridge" %in% model) {
    pt1 <- proc.time() 
    if(exists("Xtest")){ #
      ridgeres <- regress(Xtest, ytest, method = c("ridge"), p = 0)
    } else {
      ridgeres <- regress(as.matrix(X), as.matrix(y), method = c("ridge"), p = 0) }
    ridgeTime <- as.numeric(proc.time() - pt1)[3]
  }
# getting bayesian #
if ("blasso" %in% model) {
  pt1 <- proc.time()
  if(exists("Xtest")){
    blassores <- blasso(Xtest, ytest, RJ=FALSE, M=M, T=25000, verb = 0)
  } else {blassores <- blasso(as.matrix(X), as.matrix(y), RJ=FALSE, M=M, T=25000, verb = 0)}
  blasTime <- as.numeric(proc.time() - pt1)[3]
}
if ("bhs" %in% model) {
  pt1 <- proc.time()
  if(exists("Xtest")){bhsres <- bhs(Xtest, ytest, RJ=FALSE, M=M, T=25000, verb = 0)} else { 
    bhsres <- bhs(as.numeric(X), as.numeric(y), RJ=FALSE, M=M, T=25000, verb = 0)}
  bhsTime <- as.numeric(proc.time() - pt1)[3]
}
if ("bridge" %in% model) {
  pt1 <- proc.time()
  if(exists("xtest")){bridgeres <- bridge(Xtest, ytest, RJ=FALSE, M=M, T=25000, verb = 0)} else {
    bridgeres <- bridge(as.numeric(X), as.numeric(y), RJ=FALSE, M=M, T=25000, verb = 0)}
  briTime <- as.numeric(proc.time() - pt1)[3]
}
if ("susie" %in% model ) {
  pt1 <- proc.time()
  if(exists("Xtest")){
    susres <- susie(Xtest, ytest, L=M)} else { 
      susres <- susie(as.matrix(X[,-1]), as.matrix(y), L=M) }
  susTime <- as.numeric(proc.time() - pt1)[3]
}

  

#predict y from models
if (withPred == TRUE) {
  if(exists("lmres")){
    lmpred <- lmres$coefficients[1] + Xval %*% lmres$coefficients[-1] 
    #note: R appears to still not have a good way to handle predictions of different shape than initial data, so I tried to do it the old fashioned way
    lmpredRMSE <- sqrt(mean((summary(lm(lmpred~yval))$residuals)^2))
    lmpredR2 <- summary(lm(lmpred~yval))$adj.r.squared
    lmpredAIC <- AIC(lmpred)
    lmpredcor <- cor.test(lmpred, yval)
    lmPStats <- c(lmpredRMSE, lmpredR2, lmpredAIC, lmpredcor)
  }
  if(exists("olsres")){
    olspred <- olsres$b[1] + Xval %*% olsres$b[-1]
    olspredRMSE <- sqrt(mean((summary(lm(olspred~yval))$residuals)^2))
    olspredR2 <- summary(lm(olspred~yval))$adj.r.squared
    olspredAIC <- AIC(olspred)
    olspredcor <- cor.test(olspred, yval) # are lm() and regress(method = "ols") equivalent? it looks like it
    olsPStats <- c(olspredRMSE, olspredR2, olspredAIC, olspredcor)
  }
  if(exists("plsres")){
    plspred <- plsres$b[1] + Xval %*% plsres$b[-1]
    plspredRMSE <- sqrt(mean((summary(lm(plspred~yval))$residuals)^2))
    plspredR2 <- summary(lm(plspred~yval))$adj.r.squared
    plspredAIC <- AIC(plspred)
    plspredcor <- cor.test(plspred, yval) 
    plsPStats <- c(plspredRMSE, plspredR2, plspredAIC, plspredcor)
  }
  if(exists("pcrres")){
    pcrpred <- pcrres$b[1] + Xval %*% pcrres$b[-1]
    pcrpredRMSE <- sqrt(mean((summary(lm(pcrpred~yval))$residuals)^2))
    pcrpredR2 <- summary(lm(pcrpred~yval))$adj.r.squared
    pcrpredAIC <- AIC(pcrpred)
    pcrpredcor <- cor.test(pcrpred, yval)
    pcrPStats <- c(pcrpredRMSE, pcrpredR2, pcrpredAIC, pcrpredcor)
  }
  if(exists("larres")){
    larpred <- larres$b[1] + Xval %*% larres$b[-1]
    larpredRMSE <- sqrt(mean((summary(lm(larpred~yval))$residuals)^2))
    larpredR2 <- summary(lm(larpred~yval))$adj.r.squared
    larpredAIC <- AIC(larpred)
    larpredcor <- cor.test(larpred, yval) 
    larPStats <- c(larpredRMSE, larpredR2, larpredAIC, larpredcor)
  }
  if(exists("lassores")){
    lassopred <- coef(lassores)[1] + Xval %*% coef(lassores)[-1]
    lassopredRMSE <- sqrt(mean((summary(lm(lassopred~yval))$residuals)^2))
    lassopredR2 <- summary(lm(lassopred~yval))$adj.r.squared
    lassopredAIC <- AIC(lassopred)
    lassopredcor <- cor.test(lassopred, yval) 
    lassoPStats <- c(lassopredRMSE, lassopredR2, lassopredAIC, lassopredcor)
  }
  if(exists("ridgeres")){
    ridgepred <- coef(ridgeres)[1] + Xval %*% coef(ridgeres)[-1]
    ridgepredRMSE <- sqrt(mean((summary(lm(ridgepred~yval))$residuals)^2))
    ridgepredR2 <- summary(lm(ridgepred~yval))$adj.r.squared
    ridgepredAIC <- AIC(ridgepred)
    ridgepredcor <- cor.test(ridgepred, yval) 
    ridgePStats <- c(ridgepredRMSE, ridgepredR2, ridgepredAIC, ridgepredcor)
  }
  if(exists("blassores")){
    blassopred <- mean(blassores$mu[10001:25000]) + Xval[,-1] %*% colMeans(blassores$beta[10001:25000,])
    blaspredRMSE <- sqrt(mean((summary(lm(blassopred~yval))$residuals)^2))
    blaspredR2 <- summary(lm(blassopred~yval))$adj.r.squared
    blaspredAIC <- AIC(blaspred)
    blaspredcor <- cor.test(blassopred, yval) 
    blasPStats <- c(blaspredRMSE, blaspredR2, blaspredAIC, blaspredcor)
  }
  if(exists("bridgeres")){
    bridgepred <- mean(bridgeres$mu[10001:25000]) + Xval[,-1] %*% colMeans(bridgeres$beta[10001:25000,]) 
    bripredRMSE <- sqrt(mean((summary(lm(bridgepred~yval))$residuals)^2))
    bripredR2 <- summary(lm(bridgepred~yval))$adj.r.squared
    bripredAIC <- AIC(bridgepred)
    bripredcor <- cor.test(bridgepred, yval) 
    briPStats <- c(bripredRMSE, bripredR2, bripredAIC, bripredcor)
  }
  if(exists("bhsres")){
    bhspred <- mean(bhsres$mu[10001:25000]) + Xval[,-1] %*% colMeans(bhsres$beta[10001:25000,]) 
    bhspredRMSE <- sqrt(mean((summary(lm(bhspred~yval))$residuals)^2))
    bhspredR2 <- summary(lm(bhspred~yval))$adj.r.squared
    bhspredAIC <- AIC(bhspred)
    bhspredcor <- cor.test(bhspred, yval) 
    bhsPStats <- c(bhspredRMSE, bhspredR2, bhspredAIC, bhspredcor)
  }
  if(exists("susres")){
    suspred <- coef(susres)[1] + Xval[,-1] %*% coef(susres)[-1]
    #going back to my note above, predict doesn't work with a differently shaped NewData. I tried doubling my oos so it would be the same size, but it looks like pred still used x values from my original set.
    suspredRMSE <- sqrt(mean((summary(lm(suspred~yval))$residuals)^2))
    suspredR2 <- summary(lm(suspred~yval))$adj.r.squared
    suspredAIC <- AIC(suspred)
    suspredcor <- cor.test(suspred, yval) 
    susPStats <- c(suspredRMSE, suspredR2, suspredAIC, suspredcor)
  }
}
  
# it's quicker and more efficient to output performance metrics rather than model parameters, but I'll leave that as an option
  

  

#get output
out_test <- list(if(exists("ytest")){ytest}, if(exists("yval")){yval}, 
            if(exists("Xtest")){Xtest}, if(exists("Xval")){Xval}, 
            if(exists("lmres")){coef(lmres)}, if(exists("lmpred")){c(lmpred, lmpredRMSE)},
            #add ols, pls, pcr, lar
            
             if(exists("lassopred")){c(lassopred, lassoRMSE)},
            if(exists("ridgepred")){c(ridgepred, ridgeRMSE)},
            if(exists("blassopred")){c(blassopred,blassopredRMSE)}, # now, I'll have intercept and coefficients, in same form as lm and susie
             if(exists("bridgepred")){c(bridgepred, bripredRMSE)}, 
            if(exists("bhspred")){c(bhspred,bhspredRMSE)}, 
            if(exists("suspred")){c(suspred,suspredRMSE)}
)


if (pred=TRUE){
  pred_out <- 
}

if (coeff=TRUE) { # IF user wants to keep all the estimated coefficients, this is how they could do that
coeff_out <- list(if(exists("lmres")){coef(lmres)}, if(exists("olsres")){olsres$b}, if(exists("plsres")){plsres$b}, if(exists("pcrres")){pcrres$b}, if(exists("larres")){larres$b}, if(exists("lassores")){coef(lassores)}, if(exists("ridge")){coef(ridgeres)}, if(exists("blassores")){c(mean(blassores$mu[10001:25000]),colMeans(blassores$beta[10001:25000,]))}, if(exists("bridgeres")){c(mean(bridgeres$mu[10001:25000]),colMeans(bridgeres$beta[10001:25000,]))}, if(exists("bhsres")){c(mean(bhsres$mu[10001:25000]),colMeans(bhsres$beta[10001:25000,]))}, if(exists("susres")){coef(susres)})
  names(coeff_out) <- c("lmBeta", "olsBeta", "plsBeta", "pcrBeta", "larBeta", "lassoBeta", "ridgeBeta", "blasBeta", "bridgeBeta", "bhsBeta", "susBeta")
coeff_out <- coeff_out[!sapply(coeff_out,is.null)]
}
  

names(out) <- c("ytest", "yval", "Xtest", "Xval", "lmPred", "lmpredRMSE",  "lassoPred", "lassopredRMSE", "ridgePred", "ridgepredRMSE", "blasPred", "blaspredRMSE", "bridgePred", "bridgepredRMSE", "bhsPred", "bhspredRMSE", "susBeta", "susPred", "suspredRMSE")
out <- out[!sapply(out,is.null)]
return(out)  
}




###################################################################

### now run it ###
### **** will try rewriting this so that I don't need to load all my massive simulated data into memory at the same time **** ###
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
  resData <- as.data.frame(read.csv(file=resFiles[l], header = T)[,2])
  predData <- read.csv(file=predFiles[l],header=T)[,2:801]
  # check that predictor datasets match length of response datasets; if not, send warning
  if (nrow(resData) != nrow(predData)){
    print("Predictors and variables do not match. Check input.")
  }
  
  # if I want to use my already simulated out of sample observations
  if(exists("yOOSFiles")){
    yOOS <- as.data.frame(read.csv(file=yOOSFiles[l],header=T)[,2])
  }
  if(exists("XOOSFiles")){
    XOOS <- read.csv(file=XOOSFiles[l],header=T) [,2:801]
  }
       
  # fit the models, make predictions, return model diagnostics and coefficients (if requested)
  modOut <- fitMod(predData, resData, 15, yOOs, XOOS, model = c("lm", "pls", "pcr", "lar", "lasso", "ridge", "elnet", "blasso", "bridge", "bhs", "susie"), withPred=TRUE, coeff=FALSE) 
  
  # save output
  for (l in 1:length(modOut)){ 
    write.csv(stepOut[l], file = paste("Model", rep, "_", names(modOut)[l], ".csv", sep = ""))
  }
  
}



###############################
# which model did best on the excluded data?
# note: this code isn't working the way I want it to, yet
# this belongs in my subsequent script where I do comparisons
if (length(model) > 1 & withPred == TRUE){
  leastErrorPred <- c("Linear","Susie","Bayesian Lasso", "Horseshoe Lasso", "Ridge Lasso")[which.min(c(if(exists("lmpredRMSE")){lmpredRMSE}, if(exists("suspredRMSE")){suspredRMSE}, if(exists("blaspredRMSE")){blaspredRMSE}, if(exists("bhspredRMSE")){bhspredRMSE}, if(exists("bripredRMSE")){bripredRMSE}))]
}	

if (length(model) > 1){ # and which model was quickest?
  leastTime <- c("Linear","Susie","Bayesian Lasso", "Horseshoe Lasso", "Ridge Lasso")[which.min(c(if(exists("lmTime")){lmTime}, if(exists("susTime")){susTime}, if(exists("blasTime")){blasTime}, if(exists("bhsTime")){bhsTime}, if(exists("briTime")){briTime}))]
}	



############################################################3
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




