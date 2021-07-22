# Test models to see how they perform

  # things I would like to add but haven't yet:
    # BSLMM
    # RJ vs non-RJ

### load modules ###
#module load swset/2018.05  gcc/7.3.0 py-matplotlib/2.2.2-py36 r/4.0.5-py27

#try using -m 50M with salloc?

### load libraries ###
library(susieR)
library(monomvn)
library(glmnet)
library(dplyr)
library(reshape2)
library(clusterGeneration) #a more flexible way to generate covariance matrices

### get the data ###
resFiles <- c(list.files(pattern = "simy.csv", full.names=TRUE)) # create list of files containing my response variables
resData <- lapply(resFiles, function(x){read.csv(file=x,header=T)})
resAll <- Reduce(function(x,y){full_join(x,y, by = "X")}, resData) # combine all response files into a single dataframe
colnames(resAll) <- c("X", paste("rep",1:(ncol(resAll)-1), sep = "")) # a dataframe with rnows of predicted responses (individuals) over ncolumns + 1 (reps) 

  #if working with a large number of predictors, this can run out of 
predFiles <- c(list.files(pattern = "simX.csv", full.names=TRUE)) # create list of files containing my predictor variables
predData <- lapply(predFiles, function(x){read.csv(file=x,header=T)}) # a list of dataframes... values for multiple predictors for every individual across all reps (list indices)
names(predData) <- c(paste("rep",1:length(predData), sep = ""))


# check that predictor datasets match length of response datasets; if not, send warning
if (nrow(resAll) != nrow(as.data.frame(predData[1])) | ncol(resAll) -1 != length(predData)){
  print("Predictors and variables do not match. Check input.")
}


# if I want to use my already simulated out of sample observations
yOOSFiles <- c(list.files(pattern = "yOOS.csv", full.names=TRUE)) 
yOOSData <- lapply(yOOSFiles, function(x){read.csv(file=x,header=T)})
yOOSAll <- Reduce(function(x,y){full_join(x,y, by = "X")}, yOOSData)
colnames(yOOSAll) <- c("X", paste("rep",1:ncol(yOOSAll), sep = "")) 

XOOSFiles <- c(list.files(pattern = "XOOS.csv", full.names=TRUE)) 
XOOSData <- lapply(XOOSFiles, function(x){read.csv(file=x,header=T)}) 
names(XOOSData) <- c(paste("rep",1:length(XOOSData), sep = ""))



### define modeling function ###
fitMod <- function(X, y, M=NULL, yOOS=NULL, XOOS=NULL, model = c("lm", "susie", "blasso", "bridge", "bhs", "lasso",), withPred=TRUE){
  X <- X
  y <- y
  if(!is.null(M)){ M <- M }# can set max expected non-zero predictors; if sparsity is expected, lower values are faster. 

  # I saved my out of sample observations separately, so I will read that in and use
    #but I can create an test sample even without the pre-simulated set
if (withPred == TRUE){
    if (!is.null(yOOS) & !is.null(XOOS)) { # need both
    Xtest <- X
    Xval <- XOOS
    ytest <- y
    yval <- yOOS
    
  } else {   # create an out of sample 
    if (length(y)%%10 == 0){  
      valN <- length(y)/10   # I've split the sample 9:1, which works for any dataset with nObservations divisible by 10
    } else {valN <-  round(length(y)/10, 0)}   # for samples with a different length, I'll have to round 
    whichVal <- sample(1:length(y), valN)
    Xtest <- X[-(whichVal + 1),] 
    Xval <- X[(whichVal + 1),]
    ytest <- y[-(whichVal)]
    yval <- y[whichVal]
    
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
  if(exists(Xtest)){
    lmres <- lm(ytest ~ Xtest) } else {
      lmres <- lm(y ~ X) }
  lmTime <- proc.time() - pt1
}
if ("susie" %in% model ) {
  pt1 <- proc.time()
  if(exists(Xtest)){
    susres <- susie(as.matrix(Xtest[,-1]), as.matrix(ytest), L=M)} else { 
      susres <- susie(as.matrix(X[,-1]), as.matrix(y), L=M) }
  susTime <- proc.time() - pt1
}
if ("blasso" %in% model) {
  pt1 <- proc.time()
  if(exists(Xtest)){
    blassores <- blasso(Xtest[,-1], ytest, RJ=FALSE, M=M, T=25000, verb = 0)
  } else {blassores <- blasso(X[,-1], y, RJ=FALSE, M=M, T=25000, verb = 0)}
  blasTime <- proc.time() - pt1
}
if ("bhs" %in% model) {
  pt1 <- proc.time()
  if(exists(Xtest)){bhsres <- bhs(Xtest[,-1], ytest, RJ=FALSE, M=M, T=25000, verb = 0)} else { bhsres <- bhs(X[,-1], y, RJ=FALSE, M=M, T=25000, verb = 0)}
  bhsTime <- proc.time() - pt1
}
if ("bridge" %in% model) {
  pt1 <- proc.time()
  if(exists(xtest)){bridgeres <- bridge(Xtest[,-1], ytest, RJ=FALSE, M=M, T=25000, verb = 0)} else {bridgeres <- bridge(X[,-1], y, RJ=FALSE, M=M, T=25000, verb = 0)}
  briTime <- proc.time() - pt1
}
if ("lasso" %in% model){
  pt1 <- proc.time()
  if(exists(xtest)){lassores <- glmnet(Xtest[,-1], ytest, alpha=1,verb = 0)} else {lassores <- glmnet(X[,-1], y, verb = 0)}
  lasTime <- proc.time() - pt1
}
  if ("ridge" %in% model){
    pt1 <- proc.time()
    if(exists(xtest)){ridgeres <- glmnet(Xtest[,-1], ytest, alpha=0, verb = 0)} else {ridgeres <- glmnet(X[,-1], y, RJ=FALSE, M=M, T=25000, verb = 0)}
    ridgeTime <- proc.time() - pt1
  }
  if ("elnet" %in% model){
    pt1 <- proc.time()
    if(exists(xtest)){elnetres <- glmnet(Xtest[,-1], ytest, alpha=.5, verb = 0)} else {elnetres <- glmnet(X[,-1], y, verb = 0)} #I can specify the elasticnet mixing parameter between 1 and 0; I've defaulted to midway between
    elnetTime <- proc.time() - pt1
  }
#I don't fully understand glmnet's "relax" option, so I've left it out. documentation says that this feature is best if fitting the original glmnet model takes a long time.
  

#predict y from models
  #### need to add glmnet options
if (withPred == TRUE) {
  
  if(exists("lmres")){
    lmpred <- lmres$coefficients[1] + Xval[,-1] %*% lmres$coefficients[-1] 
    #note: R appears to still not have a good way to handle predictions of different shape than initial data, so I tried to do it the old fashioned way
    lmspredRMSE <- sqrt(mean((summary(lm(lmpred~yval))$residuals)^2))
  }
  if(exists("susres")){
    suspred <- coef(susres)[1] + Xval[,-1] %*% coef(susres)[-1]
    #going back to my note above, predict doesn't work with a differently shaped NewData. I tried doubling my oos so it would be the same size, but it looks like pred still used x values from my original set.
    suspredRMSE <- sqrt(mean((summary(lm(suspred~yval))$residuals)^2))
  }
  if(exists("blassores")){
    blassopred <- mean(blassores$mu[10001:25000]) + Xval[,-1] %*% colMeans(blassores$beta[10001:25000,])
    blaspredRMSE <- sqrt(mean((summary(lm(blassopred~yval))$residuals)^2))
  }
  if(exists("bhsres")){
    bhspred <- mean(bhsres$mu[10001:25000]) + Xval[,-1] %*% colMeans(bhsres$beta[10001:25000,]) 
    bhspredRMSE <- sqrt(mean((summary(lm(bhspred~yval))$residuals)^2))
  }
  if(exists("bridgeres")){
    bridgepred <- mean(bridgeres$mu[10001:25000]) + Xval[,-1] %*% colMeans(bridgeres$beta[10001:25000,]) 
    bripredRMSE <- sqrt(mean((summary(lm(bripred~yval))$residuals)^2))
  }
  if(exists("lassores")){
    lassopred <- coef(lassores)[1] + Xval[,-1] %*% coef(lassores)[-1]
    lassopredRMSE <- sqrt(mean((summary(lm(lassopred~yval))$residuals)^2))
  }
  if(exists("ridgeres")){
    ridgepred <- coef(ridgeres)[1] + Xval[,-1] %*% coef(ridgeres)[-1]
    ridgepredRMSE <- sqrt(mean((summary(lm(ridgepred~yval))$residuals)^2))
  }
  if(exists("elasres")){
    elaspred <- coef(elasres)[1] + Xval[,-1] %*% coef(elasres)[-1]
    elaspredRMSE <- sqrt(mean((summary(lm(elaspred~yval))$residuals)^2))
  }
  
}
  
# I think it would be more efficient to output performance metrics rather than model parameters, but I'm still working out how best to get comparable metrics for all the models

# which model did best on the excluded data?
  # note: this code isn't working the way I want it to, yet
if (length(model) > 1 & withPred == TRUE){
    leastErrorPred <- c("Linear","Susie","Bayesian Lasso", "Horseshoe Lasso", "Ridge Lasso")[which.min(c(if(exists("lmpredRMSE")){lmpredRMSE}, if(exists("suspredRMSE")){suspredRMSE}, if(exists("blaspredRMSE")){blaspredRMSE}, if(exists("bhspredRMSE")){bhspredRMSE}, if(exists("bripredRMSE")){bripredRMSE}))]
}	
  
if (length(model) > 1){ # and which model was quickest?
    leastTime <- c("Linear","Susie","Bayesian Lasso", "Horseshoe Lasso", "Ridge Lasso")[which.min(c(if(exists("lmTime")){lmTime}, if(exists("susTime")){susTime}, if(exists("blasTime")){blasTime}, if(exists("bhsTime")){bhsTime}, if(exists("briTime")){briTime}))]
}	
  

#get output
out <- list(if(exists("ytest")){ytest}, if(exists("yval")){yval}, 
            if(exists("Xtest")){Xtest}, if(exists("Xtest")){Xtest}, 
            if(exists("lmres")){coef(lmres)}, if(exists("lmpred")){c(lmpred, lmpredRMSE)}, 
            if(exists("susres")){coef(susres)}, if(exists("suspred")){c(suspred,suspredRMSE)},
            if(exists("blassores")){c(mean(blassores$mu[10001:25000]),colMeans(blassores$beta[10001:25000,]))}, if(exists("blassopred")){c(blassopred,blassopredRMSE)},
            # now, I'll have intercept and coefficients, in same form as lm and susie
            if(exists("bhsres")){c(mean(bhsres$mu[10001:25000]),colMeans(bhsres$beta[10001:25000,]))}, if(exists("bhspred")){c(bhspred,bhspredRMSE)},
            if(exists("bridgeres")){c(mean(bridgeres$mu[10001:25000]),colMeans(bridgeres$beta[10001:25000,]))}, if(exists("bridgepred")){c(bridgepred, bripredRMSE)}
            if(exists("lassores")){coef(lassores)}, if(exists("lassopred")){c(lassopred, lassoRMSE)},
            if(exists("ridge")){coef(ridgeres)}, if(exists("ridgepred")){c(ridgepred, ridgeRMSE)},
            if(exists("elasticnet")){coef(elasres)}, if(exists("elaspred")){c(elaspred, elasRMSE)}
)

names(out) <- c("ytest", "yval", "Xtest", "Xval", "lmBeta", "lmPred", "lmpredRMSE", "susBeta", "susPred", "suspredRMSE", "blasBeta", "blasPred", "blaspredRMSE", "bhsBeta", "bhsPred", "bhspredRMSE", "bridgeBeta", "bridgePred", "bridgepredRMSE", "lassoBeta", "lassoPred", "lassopredRMSE","ridgeBeta", "ridgePred", "ridgepredRMSE","elasBeta", "elasPred", "elaspredRMSE")
out <- out[!sapply(out,is.null)]
return(out)  
}


### now run it ###
for (rep in 1:length(resData)){ # fit models for each rep
  Xrep <- as.data.frame(predData[rep])
  yrep <- resAll[,rep]
  yOOS <- yOOSAll[,rep]
  XOOS <- as.data.frame(XOOS[rep])
  modOut <- fitMod(Xrep, yrep, 15, yOOs, XOOS, model = c("lm", "susie", "lasso", "bridge", "bhs"), withPred=TRUE)
  # save output
  for (l in 1:length(modOut)){ 
    write.csv(stepOut[l], file = paste("Model", rep, "_", names(modOut)[l], ".csv", sep = ""))
  }
}

