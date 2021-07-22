# Simulating datasets to compete models on
#module load swset/2018.05  gcc/7.3.0 py-matplotlib/2.2.2-py36 r/4.0.5-py27


# Load Libraries
library(susieR)
library(rethinking)
library(MASS) # because I don't have rethinking on the teton server
library(monomvn)
library(reshape2)
library(clusterGeneration) #a more flexible way to generate covariance matrices
set.seed(6721) # for n1000p800fewcorrvar10, seed = 6721


###for testing purposes, delete when functional###
  n <- 1000 
  p <- 800
  Var <- 10 #can change the allowable covariation!
###

  
# I want to run this 500 times to see how well each of these models do under many effect, non-sparse situations with correlated variables of variance 10  
  #this is really inefficient as is; should make a way to not calculate OOS variables if I'm not trying to compare predictions
  
  
sparseSim <- function(n, p, effect = c("many", "few"), corr = c("no", "yes"), covMeth = c("onion"), Var, sparse = FALSE){
  n <- n #observations
  p <- p #predictors
  Var <- Var

#create correlation matrix for predictors

#this method is pretty slow compared to other options but allows for more control
  #unifcorrmat is deathly slow. maybe stick to onion as the method.
Sigma <- genPositiveDefMat("onion", rangeVar = c(1,Var), dim=p)$Sigma
  #changing the rangeVar, I can change the allowable range for variances
# it's possible to change the eta parameter for generating the correlation matrix, either with the c-vine method or the onion method
# I left it with the default eta=1

  
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
# To be fair, these don't seem like equivalent comparisons. the code for the correlated dataset generates a range of response variables 3x that of the non-correlated random matrix.
  

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
  

#adjust sparsity (the above code for Beta approximates sparsity as well, this is just more explicit)
if (sparse == TRUE) { # if most predictors really don't have any effect, they should == 0
  beta <- beta
  beta[sample(c(1:p), p*(9/10), replace = F)] <- 0 # here, only 1/10 of predictors are allowed to be non-zero
}


#simulate observations
yfull <- Xfull %*% beta + rnorm(3*n/2) 
y <- yfull[1:n,] # this is my smaller, "sampled" population



#get output
out <- list(y, X, beta, 
            if(exists("yfull")){yfull[(n+1):nrow(yfull),]}, 
            if(exists("Xfull")){Xfull[(n+1):nrow(Xfull),]}
)

names(out) <- c("simy", "simX", "simBeta", "yOOS", "XOOS")
out <- out[!sapply(out,is.null)]
return(out)  
}


### LOOP CODE ###

AllOutput <- list()
#for (obs in c(100, 500, 1000, 5000)){ 
  # I still need to write this in a way that puts the output in different folders if I were to loop across different numbers of samples or variables
for (i in 1:500) {
  stepOut <- sparseSim(1000, 800, effect = "few", corr = "yes", Var = 10)
  repname <- paste("rep", i, sep = "")
  print(repname)
  for (l in 1:length(stepOut)){
    write.csv(stepOut[l], file = paste(repname, "_", names(stepOut)[l], ".csv", sep = ""))
  }
}
#}
#save(AllOutput, file = "filename.Rdata")
    # added a write step within the loop in case I hit a wall limit
    # writing reps individually also allows me to stay under the default memory limit
