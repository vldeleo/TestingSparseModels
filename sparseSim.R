# Competing sparse model approaches across different situations

# Load Libraries
library(susieR)
library(rethinking)
library(MASS) # because I don't have rethinking on the teton server
library(monomvn)
library(reshape2)
library(ggplot2)
library(clusterGeneration) #a more flexible way to generate covariance matrices
library(ggpubr) # for plotting in multiple panes
library(qgraph) # because I'd like cool viz of my correlated predictors
set.seed(518) #used this for first set
              # second set (with many effects) I set.seed(527) and then after 30 set.seed(528)

# Simulate Data
  #I think I'll want this to be a function that I can choose options
  #correlation: high (0-1), low (-.4-.4), no (0)
    #ncorr: many (50% predictors), few (20% predictors)
  #effect: high (), mostlylow (see beta definition below), low ()
    #neffect: many (1/5), few (1/50) (question of sparsity)
  #how many observations?
  #how many covariates?

###for testing purposes, delete when functional###
  n <- 1000 
  p <- 800
  Var <- 10 #can change the allowable covariation!
###

  
  
# I want to run this 500 times to see how well each of these models do under many effect, non-sparse situations with correlated variables of variance 10  
  #this is really inefficient as is; should make a way to not calculate OOS variables if I'm not trying to compare predictions
  

sparseSim <- function(n, p, effect = c("many", "few"), corr = c("no", "yes"), covMeth = c("onion"), Var, sparse = FALSE, model = c("lm", "susie", "blasso", "bhs", "bridge"), M, withPred = TRUE){
  n <- n #observations
  p <- p #predictors
  Var <- Var
  M <- M

#create correlation matrix for predictors

#this method is pretty slow compared to other options but allows for more control
  #unifcorrmat is deathly slow. maybe stick to onion as the method.
Sigma <- genPositiveDefMat("onion", rangeVar = c(0,Var), dim=p)$Sigma
  #changing the rangeVar, I can change the allowable range for variances
  
#### Originally I wanted an option to specify how /many/ predictors are correlated as well as how /strongly/ predictors are correlated
  # I still haven't managed that
  

#simulate values of X
#to test predictive ability, let me simulate an OOS population
  
  if (corr == "no"){
    Xfull <- matrix(rnorm((3*n/2)*p, 0, 1),nrow=(3*n/2),ncol=p)
  } else if (corr == "yes") {   # Sigma will be different for high/low based on code above
    mu <- rnorm(p, 0, 1)
    Xfull <- MASS::mvrnorm((3*n/2), mu = mu, Sigma = Sigma)  # draw predictors from multivariate normal distribution  
  } else {
  print("Invalid option given for correlation.")
  }
# To be fair, these don't seem like equivalent comparisons. the code for the correlated dataset generates a range 3x that of the non-correlated random matrix.
  

# now, go back to the appropriately sized X
  #X is dim n*p

X <- Xfull[1:n,]



  
#simulate effects
  #beta is length of p
if (effect == "many"){
  beta <- rgamma(p, 0.3, 0.25) * sample(c(-1, 1), p, replace = T)
} else if (effect == "few") {
  beta <- rgamma(p, 0.02, 0.1) * sample(c(-1, 1), p, replace = T)
}
  

#adjust sparsity (the above code for Beta causes sparsity as well, this is just more explicit)
if (sparse == TRUE) { # if most predictors really don't have any effect, they should == 0
  beta <- beta
  beta[sample(c(1:p), p*(9/10), replace = F)] <- 0 # here, only 1/10 of predictors are allowed to be non-zero
}


#simulate observations
yfull <- Xfull %*% beta + rnorm(3*n/2) 
y <- yfull[1:n,] # this is my smaller, "sampled" population


#fit model
  # I think I would like to have the option to test multiple models on the same simulated data at once
if ("lm" %in% model) {
lmres <- lm(y ~ X)
}
if ("susie" %in% model ) {
susres <- susie(X, y, L=M)
}
if ("blasso" %in% model) {
blassores <- blasso(X, y, RJ=TRUE, M=M, T=15000, verb = 0)
}
if ("bhs" %in% model) {
bhsres <- bhs(X, y, RJ=TRUE, M=M, T=15000, verb = 0)
}
if ("bridge" %in% model) {
bridgeres <- bridge(X, y, RJ=TRUE, M=M, T=15000, verb = 0)
}


#predict y from models
if (withPred == TRUE) {
if(exists("lmres")){
  lmpred <- lmres$coefficients[1] + Xfull[(n+1):nrow(Xfull),] %*% lmres$coefficients[-1] 
    #note: R appears to still not have a good way to handle predictions of different shape than initial data, so I tried to do it the old fashioned way
}
if(exists("susres")){
  suspred <- coef(susres)[1] + Xfull[(n+1):nrow(Xfull),] %*% coef(susres)[-1]
  #suspred <- predict(susres, newdata = Xfull[(n+1):nrow(Xfull),])
    #going back to my note above, predict doesn't work with a differently shaped NewData. I tried doubling my oos so it would be the same size, but it looks like pred still used x values from my original set. GAH.
}
if(exists("blassores")){
  blassopred <- mean(blassores$mu[5001:15000]) + Xfull[(n+1):nrow(Xfull),] %*% colMeans(blassores$beta[5001:15000,]) 
}
if(exists("bhsres")){
  bhspred <- mean(bhsres$mu[5001:15000]) + Xfull[(n+1):nrow(Xfull),] %*% colMeans(bhsres$beta[5001:15000,]) 
}
if(exists("bridgeres")){
  bridgepred <- mean(bridgeres$mu[5001:15000]) + Xfull[(n+1):nrow(Xfull),] %*% colMeans(bridgeres$beta[5001:15000,]) 
}
}

#get output
out <- list(y, X, beta, 
            if(exists("yfull")){yfull[(n+1):nrow(yfull)]}, if(exists("Xfull")){Xfull[(n+1):nrow(Xfull),]},
            if(exists("lmres")){coef(lmres)}, if(exists("lmpred")){lmpred},
            if(exists("susres")){coef(susres)}, if(exists("suspred")){suspred},
            if(exists("blassores")){c(mean(blassores$mu[5001:15000]),colMeans(blassores$beta[5001:15000,]))}, if(exists("blassopred")){blassopred},
            # now, I'll have intercept and coefficients, in same form as lm and susie
            if(exists("bhsres")){c(mean(bhsres$mu[5001:15000]),colMeans(bhsres$beta[5001:15000,]))}, if(exists("bhspred")){bhspred},
            if(exists("bridgeres")){c(mean(bridgeres$mu[5001:15000]),colMeans(bridgeres$beta[5001:15000,]))}, if(exists("bridgepred")){bridgepred}
            )

names(out) <- c("simy", "simX", "simBeta", "yOOS", "XOOS", "lmBeta", "lmPred", "susBeta", "susPred", "blasBeta", "blasPred", "bhsBeta", "bhsPred", "bridgeBeta", "bridgePred")
out <- out[!sapply(out,is.null)]
return(out)  
}


### LOOP CODE ###

AllOutput <- list()
# use set.seed 529 for 61-75
# set.seed(530) for 76-
for (i in 76:100) {
  stepOut <- sparseSim(500, 400, effect = "many", corr = "yes", Var = 10, model = c("lm", "susie", "blasso", "bhs", "bridge"), M = 15, withPred = TRUE)
  repname <- paste("rep", i, sep = "")
  #AllOutput[[repname]] <- stepOut
  print(repname)
  for (l in 1:length(stepOut)){
    write.csv(stepOut[l], file = paste(repname, "_", names(stepOut)[l], ".csv", sep = ""))
  }
}
save(AllOutput, file = "n1000p800fewcorrvar10m15_2021.5.21.Rdata")
# this is crazy slow and will never finish in time ):
    # added a write step within the loop in case I hit a wall limit
    # writing reps individually allows me to stay under the default memory limit for my interactive session
#for n = 1000 and p = 800, this takes about 5 min/rep
  # can do 40 reps in <4 hours, which seems like the best choice for now



#I am afraid something weird is going on with the lmres coefficients returned in my function
  #but admittedly lm() could just be a terrible way to model the data


