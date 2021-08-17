# Goal: Test sparse models against each other

  # things I would like to add but haven't yet:
    # BSLMM
    # RJ vs non-RJ
    # use the method to specify burnin instead of writing it out by hand
    

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
  # in this version, I do prediction and coefficient step right after model fitting, then save and remove from working data to save memory
  # if coefficients are saved, can later compare to true coefficients for the simulated data
  
if ("lm" %in% model) {
  pt1 <- proc.time() # I'd like to compare how long each model takes
  if(exists("Xtest")){
    lmres <- lm(ytest ~ Xtest) } else {
      lmres <- lm(as.matrix(y) ~ as.matrix(X)) }
  lmTime <- as.numeric(proc.time() - pt1)[3]
    # since I do have the /real/ coefficients, I could test my model accuracy here
  lmRMSE <- sqrt(mean((resid(lmres))^2))
  lmMAE <- mean(abs(resid(lmres)))
  lmAIC <- AIC(lmres)
  lmStats <- c(lmTime, lmRMSE, lmMAE, lmAIC)
  write.csv(lmStats, file = paste("ModelResults/Rep", l, "_lmModStats.csv", sep = ""))
  #predict y from models
  if (withPred==TRUE) {
    lmpred <- lmres$coefficients[1] + Xval %*% lmres$coefficients[-1] 
    #note: R appears to not have a good way to handle predictions of different shape than initial data
    lmpredRMSE <- sqrt(mean((resid(lm(lmpred~yval)))^2))
    lmpredR2 <- summary(lm(lmpred~yval))$adj.r.squared
    lmpredcor <- cor.test(lmpred, yval)
    lmPStats <- c(lmpredRMSE, lmpredR2, lmpredcor)
    write.csv(lmPStats, file = paste("ModelResults/Rep", l, "_lmPredStats.csv", sep = ""))
    rm(lmpred, lmpredRMSE, lmpredR2, lmpredcor, lmPStats)
  }
  if (coeff==TRUE) {
    write.csv(coef(lmres), file = paste("ModelResults/Rep", l, "_lmModCoeff.csv", sep = ""))
  }
  rm(lmres, lmTime, lmRMSE, lmMAE, lmAIC, lmStats)
  gc()
}
if ("ols" %in% model) {
  pt1 <- proc.time() 
  if(exists("Xtest")){ # p = 0 forces regress to use my chosen method even when there are more columns than rows. can choose a value 0-1 to specify when to use lsr
    olsres <- regress(Xtest, ytest, method = c("lsr"), p = 0) 
    } else {
      olsres <- regress(as.matrix(X), as.matrix(y), method = c("lsr"), p = 0) }
  olsTime <- as.numeric(proc.time() - pt1)[3]
  if(exists("Xtest")){olsResid <- c(olsres$b[1] + as.matrix(Xtest) %*% olsres$b[-1]) - ytest} else {
    olsResid <- mean(olsres$b[1] + as.matrix(X) %*% olsres$b[-1]) - y}
  olsRMSE <- sqrt(mean((olsResid)^2))
  olsMAE <- mean(abs(olsResid))
  olsAIC <- AIC(olsres)
  olsStats <- c(olsTime, olsRMSE, olsMAE, olsAIC)
  write.csv(olsStats, file = paste("ModelResults/Rep", l, "_olsModStats.csv", sep = ""))
    # possibility here to compare model coef to /real/ coeff 
  if (withPred==TRUE){
    olspred <- olsres$b[1] + Xval %*% olsres$b[-1]
    olspredRMSE <- sqrt(mean((summary(lm(olspred~yval))$residuals)^2))
    olspredR2 <- summary(lm(olspred~yval))$adj.r.squared
    olspredcor <- cor.test(olspred, yval) # are lm() and regress(method = "ols") equivalent? it looks like it
    olsPStats <- c(olspredRMSE, olspredR2, olspredcor)
    write.csv(olsPStats, file = paste("ModelResults/Rep", l, "_olsPredStats.csv", sep = ""))
    rm(olspred, olspredRMSE, olspredR2, olspredcor, olsPStats)
  }
  if (coeff==TRUE){ 
    write.csv(olsres$b, file = paste("ModelResults/Rep", l, "_olsModCoeff.csv", sep = ""))
  }
  rm(olsres, olsTime, olsRMSE, olsMAE, olsAIC, olsStats)
  gc()
}
if ("pls" %in% model) {
  pt1 <- proc.time() 
  if(exists("Xtest")){ # can set max number of pc to consider
    plsres <- regress(Xtest, ytest, method = c("plsr"), p = 0, ncomp.max = 10)
    } else {
      plsres <- regress(as.matrix(X), as.matrix(y), method = c("plsr"), p = 0, ncomp.max = 10) }
  plsTime <- as.numeric(proc.time() - pt1)[3]
  if(exists("Xtest")){plsResid <- c(plsres$b[1] + as.matrix(Xtest) %*% plsres$b[-1]) - ytest} else {
    plsResid <- mean(plsres$b[1] + as.matrix(X) %*% plsres$b[-1]) - y}
  plsRMSE <- sqrt(mean((plsResid)^2))
  plsMAE <- mean(abs(plsResid))
  plsAIC <- NA #AIC(plsres)
  plsStats <- c(plsTime, plsRMSE, plsMAE, plsAIC)
  write.csv(plsStats, file = paste("ModelResults/Rep", l, "_plsModStats.csv", sep = ""))
  if (withPred==TRUE) {
    plspred <- plsres$b[1] + Xval %*% plsres$b[-1]
    plspredRMSE <- sqrt(mean((summary(lm(plspred~yval))$residuals)^2))
    plspredR2 <- summary(lm(plspred~yval))$adj.r.squared
    plspredcor <- cor.test(plspred, yval) 
    plsPStats <- c(plspredRMSE, plspredR2, plspredcor)
    write.csv(plsPStats, file = paste("ModelResults/Rep", l, "_plsPredStats.csv", sep = ""))
    rm(plspred, plspredRMSE, plspredR2, plspredcor, plsPStats)
  }
  if (coeff==TRUE) {
    write.csv(plsres$b, file = paste("ModelResults/Rep", l, "_plsModCoeff.csv", sep = ""))
  }
  rm(plsres, plsTime, plsRMSE, plsMAE, plsAIC, plsStats)
  gc()
}
if ("pcr" %in% model) {
  pt1 <- proc.time() 
  if(exists("Xtest")){ 
    pcrres <- regress(Xtest, ytest, method = c("pcr"), p = 0, ncomp.max = 10)
  } else {
    pcrres <- regress(as.matrix(X), as.matrix(y), method = c("pcr"), p = 0, ncomp.max = 10) }
    pcrTime <- as.numeric(proc.time() - pt1)[3]
    if(exists("Xtest")){pcrResid <- c(pcrres$b[1] + as.matrix(Xtest) %*% pcrres$b[-1]) - ytest} else {
      pcrResid <- mean(pcrres$b[1] + as.matrix(X) %*% pcrres$b[-1]) - y}
    pcrRMSE <- sqrt(mean((pcrResid)^2))
    pcrMAE <- mean(abs(pcrResid))
    pcrAIC <- NA #AIC(pcrres)
    pcrStats <- c(pcrTime, pcrRMSE, pcrMAE, pcrAIC)
    write.csv(pcrStats, file = paste("ModelResults/Rep", l, "_pcrModStats.csv", sep = ""))
    if (withPred==TRUE) {
      pcrpred <- pcrres$b[1] + Xval %*% pcrres$b[-1]
      pcrpredRMSE <- sqrt(mean((summary(lm(pcrpred~yval))$residuals)^2))
      pcrpredR2 <- summary(lm(pcrpred~yval))$adj.r.squared
      pcrpredcor <- cor.test(pcrpred, yval)
      pcrPStats <- c(pcrpredRMSE, pcrpredR2, pcrpredcor)
      write.csv(pcrPStats, file = paste("ModelResults/Rep", l, "_pcrPredStats.csv", sep = ""))
      rm(pcrpred, pcrpredRMSE, pcrpredR2, pcrpredcor, pcrPStats)
    }
    if (coeff==TRUE) {
      write.csv(pcrres$b, file = paste("ModelResults/Rep", l, "_pcrModCoeff.csv", sep = ""))
    }
    rm(pcrres, pcrTime, pcrRMSE, pcrMAE, pcrAIC, pcrStats)
    gc()
}
if ("lar" %in% model) {
  pt1 <- proc.time() 
    if(exists("Xtest")){ 
    larres <- regress(Xtest, ytest, method = c("lar"), p = 0)
  } else {
    larres <- regress(as.matrix(X), as.matrix(y), method = c("lar"), p = 0) }
    larTime <- as.numeric(proc.time() - pt1)[3]
    if(exists("Xtest")){larResid <- c(larres$b[1] + as.matrix(Xtest) %*% larres$b[-1]) - ytest} else {
      larResid <- mean(larres$b[1] + as.matrix(X) %*% larres$b[-1]) - y}
    larRMSE <- sqrt(mean((larResid)^2))
    larMAE <- mean(abs(larResid))
    larAIC <- NA #AIC(larres)
    larStats <- c(larTime, larRMSE, larMAE, larAIC)
    write.csv(larStats, file = paste("ModelResults/Rep", l, "_larModStats.csv", sep = ""))
    if (withPred==TRUE) {
      larpred <- larres$b[1] + Xval %*% larres$b[-1]
      larpredRMSE <- sqrt(mean((summary(lm(larpred~yval))$residuals)^2))
      larpredR2 <- summary(lm(larpred~yval))$adj.r.squared
      larpredcor <- cor.test(larpred, yval) 
      larPStats <- c(larpredRMSE, larpredR2, larpredcor)
      write.csv(larPStats, file = paste("ModelResults/Rep", l, "_larPredStats.csv", sep = ""))
      rm(larpred, larpredRMSE, larpredR2, larpredcor, larPStats)
    }
    if (coeff==TRUE) {
      write.csv(larres$b, file = paste("ModelResults/Rep", l, "_larModCoeff.csv", sep = ""))
    }
    rm(larres, larTime, larRMSE, larMAE, larAIC, larStats)
    gc()   
}
if ("lasso" %in% model) {
    pt1 <- proc.time() 
    if(exists("Xtest")){ 
      lasres <- regress(Xtest, ytest, method = c("lasso"), p = 0)
    } else {
      lasres <- regress(as.matrix(X), as.matrix(y), method = c("lasso"), p = 0) }
    lasTime <- as.numeric(proc.time() - pt1)[3]
    if(exists("Xtest")){lasResid <- c(lasres$b[1] + as.matrix(Xtest) %*% lasres$b[-1]) - ytest} else {
      lasResid <- mean(lasres$b[1] + as.matrix(X) %*% lasres$b[-1]) - y}
    lasRMSE <- sqrt(mean((lasResid)^2))
    lasMAE <- mean(abs(lasResid))
    lasAIC <- NA #AIC(lasres) I'm going to have to find a method for this
    lasStats <- c(lasTime, lasRMSE, lasMAE, lasAIC)
    write.csv(lasStats, file = paste("ModelResults/Rep", l, "_lassoModStats.csv", sep = ""))
    if (withPred==TRUE) {
      lassopred <- lasres$b[1] + Xval %*% lasres$b[-1]
      lassopredRMSE <- sqrt(mean((summary(lm(lassopred~yval))$residuals)^2))
      lassopredR2 <- summary(lm(lassopred~yval))$adj.r.squared
      lassopredcor <- cor.test(lassopred, yval) 
      lassoPStats <- c(lassopredRMSE, lassopredR2, lassopredcor)
      write.csv(lassoPStats, file = paste("ModelResults/Rep", l, "_lassoPredStats.csv", sep = ""))
      rm(lassopred, lassopredRMSE, lassopredR2, lassopredcor, lassoPStats)
    }
    if (coeff==TRUE) {
      write.csv(coef(lasres), file = paste("ModelResults/Rep", l, "_lassoModCoeff.csv", sep = ""))
    }
    rm(lasres, lasTime, lasRMSE, lasMAE, lasAIC, lasStats)
    gc()
  }
if ("ridge" %in% model) {
    pt1 <- proc.time() 
    if(exists("Xtest")){ 
      ridgeres <- regress(Xtest, ytest, method = c("ridge"), p = 0)
    } else {ridgeres <- regress(as.matrix(X), as.matrix(y), method = c("ridge"), p = 0) }
    ridgeTime <- as.numeric(proc.time() - pt1)[3]
    if(exists("Xtest")){ridgeResid <- c(ridgeres$b[1] + as.matrix(Xtest) %*% ridgeres$b[-1]) - ytest} else {
      ridgeResid <- mean(ridgeres$b[1] + as.matrix(X) %*% ridgeres$b[-1]) - y}
    ridgeRMSE <- sqrt(mean((ridgeResid)^2))
    ridgeMAE <- mean(abs(ridgeResid))
    ridgeAIC <- NA #AIC(ridgeres)
    ridgeStats <- c(ridgeTime, ridgeRMSE, ridgeMAE, ridgeAIC)
    write.csv(ridgeStats, file = paste("ModelResults/Rep", l, "_ridgeModStats.csv", sep = ""))
    if (withPred==TRUE) {
      ridgepred <- ridgeres$b[1] + Xval %*% ridgeres$b[-1]
      ridgepredRMSE <- sqrt(mean((summary(lm(ridgepred~yval))$residuals)^2))
      ridgepredR2 <- summary(lm(ridgepred~yval))$adj.r.squared
      ridgepredcor <- cor.test(ridgepred, yval) 
      ridgePStats <- c(ridgepredRMSE, ridgepredR2, ridgepredcor)
      write.csv(ridgePStats, file = paste("ModelResults/Rep", l, "_ridgePredStats.csv", sep = ""))
      rm(ridgepred, ridgepredRMSE, ridgepredR2, ridgepredcor, ridgePStats)
    }
    if (coeff==TRUE) {
      write.csv(coeff(ridgeres), file = paste("ModelResults/Rep", l, "_ridgeModCoeff.csv", sep = ""))
    }
    rm(ridgeres, ridgeTime, ridgeRMSE, ridgeMAE, ridgeAIC, ridgeStats)
    gc()
  }

  # bayesian options #
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
  write.csv(blasStats, file = paste("ModelResults/Rep", l, "_blasModStats.csv", sep = ""))
  if (withPred==TRUE) {
    blassopred <- mean(blassores$mu[10001:25000]) + Xval %*% colMeans(blassores$beta[10001:25000,])
    blaspredRMSE <- sqrt(mean((resid(lm(blassopred~yval)))^2))
    blaspredR2 <- summary(lm(blassopred~yval))$adj.r.squared
    blaspredcor <- cor.test(blassopred, yval) 
    blasPStats <- c(blaspredRMSE, blaspredR2, blaspredcor)
    write.csv(blasPStats, file = paste("ModelResults/Rep", l, "_blasPredStats.csv", sep = ""))
    rm(blassopred, blaspredRMSE, blaspredR2, blaspredcor, blasPStats)
  }
  if (coeff==TRUE) {
    write.csv(c(mean(blassores$mu[10001:25000]),colMeans(blassores$beta[10001:25000,])), file = paste("ModelResults/Rep", l, "_blasModCoeff.csv", sep = ""))
  }
  rm(blassores, blasTime, blasRMSE, blasMAE, blasAIC, blasStats)
  gc()
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
  bhsAIC <- NA #AIC(bhsres) #need to find another approach here. is BIC good alternative?
  bhsStats <- c(bhsTime, bhsRMSE, bhsMAE, bhsAIC)
  write.csv(bhsStats, file = paste("ModelResults/Rep", l, "_bhsModStats.csv", sep = ""))
  if (withPred==TRUE) {
    bhspred <- mean(bhsres$mu[10001:25000]) + Xval %*% colMeans(bhsres$beta[10001:25000,]) 
    bhspredRMSE <- sqrt(mean((summary(lm(bhspred~yval))$residuals)^2))
    bhspredR2 <- summary(lm(bhspred~yval))$adj.r.squared
    bhspredcor <- cor.test(bhspred, yval) 
    bhsPStats <- c(bhspredRMSE, bhspredR2, bhspredcor)
    write.csv(bhsPStats, file = paste("ModelResults/Rep", l, "_bhsPredStats.csv", sep = ""))
    rm(bhspred, bhspredRMSE, bhspredR2, bhspredcor, bhsPStats)
  }
  if (coeff==TRUE) {
    write.csv(c(mean(bhsres$mu[10001:25000]),colMeans(bhsres$beta[10001:25000,])), file = paste("ModelResults/Rep", l, "_bhsModCoeff.csv", sep = ""))
  }
  rm(bhsres, bhsTime, bhsRMSE, bhsMAE, bhsAIC, bhsStats)
  gc()
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
  write.csv(briStats, file = paste("ModelResults/Rep", l, "_briModStats.csv", sep = ""))
  if (withPred==TRUE) {
    bridgepred <- mean(bridgeres$mu[10001:25000]) + Xval %*% colMeans(bridgeres$beta[10001:25000,]) 
    bripredRMSE <- sqrt(mean((summary(lm(bridgepred~yval))$residuals)^2))
    bripredR2 <- summary(lm(bridgepred~yval))$adj.r.squared
    bripredcor <- cor.test(bridgepred, yval) 
    briPStats <- c(bripredRMSE, bripredR2, bripredcor)
    write.csv(briPStats, file = paste("ModelResults/Rep", l, "_briPredStats.csv", sep = ""))
    rm(bridgepred, bripredRMSE, bripredR2, bripredcor, briPStats)
  }
  if (coeff==TRUE) {
    write.csv(c(mean(bridgeres$mu[10001:25000]),colMeans(bridgeres$beta[10001:25000,])), file = paste("ModelResults/Rep", l, "_bridgeModCoeff.csv", sep = ""))
  }
  rm(bridgeres, briTime, briRMSE, briMAE, briAIC, briStats)
  gc()
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
  write.csv(susStats, file = paste("ModelResults/Rep", l, "_susModStats.csv", sep = ""))
  if (withPred==TRUE) {
    suspred <- coef(susres)[1] + Xval %*% coef(susres)[-1]
    suspredRMSE <- sqrt(mean((summary(lm(suspred~yval))$residuals)^2))
    suspredR2 <- summary(lm(suspred~yval))$adj.r.squared
    suspredcor <- cor.test(suspred, yval) 
    susPStats <- c(suspredRMSE, suspredR2, suspredcor)
    write.csv(susPStats, file = paste("ModelResults/Rep", l, "_susPredStats.csv", sep = ""))
    rm(suspred, suspredRMSE, suspredR2, suspredcor, susPStats)
  }
  if (coeff==TRUE) {
    write.csv(coef(susres), file = paste("ModelResults/Rep", l, "_susModCoeff.csv", sep = ""))
  }
  rm(susres, susTime, susRMSE, susMAE, susAIC, susStats)
  gc()
}


#if creating OOS set for testing rather than using existing OOS, save the input
if (is.null(yOOS) | is.null(XOOS)){
  write.csv(c("yval", "Xval"), file = paste("ModelResults/Rep", l, "_FitInput.csv", sep = ""))
  write.csv(c("ytest", "Xtest"), file = paste("ModelResults/Rep", l, "_TestInput.csv", sep = "")) 
}
return("output saved to file")
}




###################################################################
### now run it ###
### rewrote this so that I don't load all my simulated data into memory at the same time and crash my session ###
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
  fitMod(predData, resData, 15, yOOS, XOOS, model = c("lm", "pls", "pcr", "lar", "lasso", "ridge", "blasso", "bridge", "bhs", "susie"), withPred=TRUE, coeff=FALSE) 
  #rm(list=ls()) looking for a way to clean up intermediary files so my process doesn't get killed... this might be too much, though, since it deletes my function
  gc() #this may not be enough. 
  print(l)
}

#################################################
