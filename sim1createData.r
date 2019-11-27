#this file creates the simulation datasets for the first set of simulations in the paper

library(MASS)

corr <- matrix(c(0.4,0.2,0.2,0.4), nrow=2)

#sample size
n <- 500

#it creates 1000 batches of 10 datasets, with a view to being analysed on a high performance
#computing cluster, where each instance will analyse a batch of 10 simulated datasets
nBatches <- 1000
nPerBatch <- 10

for (batch in 1:nBatches) {

  originalDataList <- list()

  for (sim in 1:nPerBatch) {

  #generate no deviation data
  z <- c(rep(0,n/2), rep(1,n/2))
  y <- mvrnorm(n=n, mu=c(2,2), Sigma=corr)
  y[z==1,2] <- y[z==1,2]+0.2
  
  #generate deviation indicator
  d <- 1*(runif(n)<0.5)
  
  #now make postdeviation data missing
  y[(d==1),2] <- NA
  
  originalDataList[[sim]] <- data.frame(z=z,y=y)

  }
  
  save(originalDataList, file=paste("./datasets/inputData_", batch, ".RData", sep=""))
}

