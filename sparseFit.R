# Test models to see how they perform

  # things I would like to add but haven't yet:
    # BSLMM
    # RJ vs non-RJ
    # non-bayesian lasso

### load libraries ###
library(susieR)
library(monomvn)
library(dplyr)
library(reshape2)
library(clusterGeneration) #a more flexible way to generate covariance matrices

### get the data ###
resFiles <- c(list.files(pattern = "simy.csv", full.names=TRUE)) # create list of files containing my response variables
resData <- lapply(resFiles, function(x){read.csv(file=x,header=T)})
resAll <- Reduce(function(x,y){full_join(x,y, by = "X")}, resData) # combine all response files into a single dataframe
colnames(resAll) <- c("X", paste("rep",1:ncol(resAll), sep = "")) # a dataframe with rnows of predicted responses (individuals) over ncolumns + 1 (reps) 

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
fitMod <- function(X, y, M=NULL, yOOS=NULL, XOOS=NULL, model = c("lm", "susie", "blasso", "bridge", "bhs"), withPred=TRUE){
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
  if(exists(Xtest)){lmres <- lm(ytest ~ Xtest) }
  else {lmres <- lm(y ~ X) }
}
if ("susie" %in% model ) {
  if(exists(Xtest)){susres <- susie(as.matrix(Xtest[,-1]), as.matrix(ytest), L=M)}
  else { susres <- susie(as.matrix(X[,-1]), as.matrix(y), L=M) }
}
if ("blasso" %in% model) {
  if(exists(Xtest)){blassores <- blasso(Xtest[,-1], ytest, RJ=FALSE, M=M, T=25000, verb = 0)} else {blassores <- blasso(X[,-1], y, RJ=FALSE, M=M, T=25000, verb = 0)}
}
if ("bhs" %in% model) {
  if(exists(Xtest)){bhsres <- bhs(Xtest[,-1], ytest, RJ=FALSE, M=M, T=25000, verb = 0)} else { bhsres <- bhs(X[,-1], y, RJ=FALSE, M=M, T=25000, verb = 0)}
}
if ("bridge" %in% model) {
  if(exists(xtest)){bridgeres <- bridge(Xtest[,-1], ytest, RJ=FALSE, M=M, T=25000, verb = 0)} else {bridgeres <- bridge(X[,-1], y, RJ=FALSE, M=M, T=25000, verb = 0)}
}


#predict y from models
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
}
  
# I think it would be more efficient to output performance metrics rather than model parameters, but I'm still working out how best to get comparable metrics for all the models

# which model did best on the excluded data?
if (length(model) > 1 & withPred == TRUE){
    leastErrorPred <- c("Linear","Susie","Bayesian Lasso")[which.min(c(if(exists("lmpredRMSE")){lmpredRMSE}, if(exists("suspredRMSE")){suspredRMSE}, if(exists("blaspredRMSE")){blaspredRMSE}))]
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
)

names(out) <- c("ytest", "yval", "Xtest", "Xval", "lmBeta", "lmPred", "lmpredRMSE", "susBeta", "susPred", "suspredRMSE", "blasBeta", "blasPred", "blaspredRMSE", "bhsBeta", "bhsPred", "bhspredRMSE", "bridgeBeta", "bridgePred", "brepredRMSE")
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

