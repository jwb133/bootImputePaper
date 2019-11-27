#the file loads up the datasets created as a result of analysing the simulations
#summary statistics are calculated and put into LaTeX tables for the paper

library(abind)
library(xtable)

load("./results/simRes_1.RData")

j2rResultsLoad <- j2rResults
marResultsLoad <- marResults

nBatches <- 1000
for (i in 2:nBatches) {
  load(paste("simRes_",i,".RData",sep=""))
  j2rResultsLoad <- abind(j2rResultsLoad,j2rResults,along=1)
  marResultsLoad <- abind(marResultsLoad,marResults,along=1)
}
rm(j2rResults,marResults)

colMeans(j2rResultsLoad[,,1])
apply(j2rResultsLoad[,,1],2,var)
colMeans(j2rResultsLoad[,,2])

numMethods <- dim(j2rResultsLoad)[2]
methodNames <- c("MI Rubin", "MI boot Rubin", "MI boot pooled percentile", "miBootPooledNormal", "Boot MI percentile", "Boot MI (M=1) percentile",
                 "bootMInormal", "bootMI1normal", "vonHippelAll", "von Hippel")

myTable <- function(results, true, methodsToPrint, ciOnly=FALSE) {
  resultsToPrint <- results[,methodsToPrint,]
  if (ciOnly==FALSE) {
    res <- array(0, dim=c(length(methodsToPrint),4))
  } else {
    res <- array(0, dim=c(length(methodsToPrint),2))
  }
  colNum <- 1
  if (ciOnly==FALSE) {
    res[,1] <- apply(resultsToPrint[,,1],2,sd)
    res[,2] <- colMeans(resultsToPrint[,,2]^0.5)
    colNum <- 3
  }
  for (i in 1:length(methodsToPrint)) {
    res[i,colNum] <- median(resultsToPrint[,i,4]-resultsToPrint[,i,3])
    res[i,colNum+1] <- 100*mean((resultsToPrint[,i,3]<true) & (resultsToPrint[,i,4]>true))
  }
  rownames(res) <- methodNames[methodsToPrint]
  if (ciOnly==FALSE) {
    colnames(res) <- c("Emp. SE", "Est. SE", "Med. CI width", "CI coverage")
  } else {
    colnames(res) <- c("Med. CI width", "CI coverage")
  }
  res
}

#first look at all methods
myTable(marResultsLoad,0.2,1:10)
myTable(j2rResultsLoad,0.1,1:10)

#now just those we will print in the paper
myTable(marResultsLoad,0.2,c(1,2,3,5,6,10))
myTable(j2rResultsLoad,0.1,c(1,2,3,5,6,10))

myTable(j2rResultsLoad,0.1,c(1,2,3,5,6,10), ciOnly=TRUE)
xtable(myTable(j2rResultsLoad,0.1,c(1,2,3,5,6,10), ciOnly=TRUE), digits=c(3,3,2))

#now print MAR and J2R side by side
cbind(myTable(marResultsLoad,0.2,c(1,2,3,5,6,10), ciOnly=TRUE), myTable(j2rResultsLoad,0.1,c(1,2,3,5,6,10), ciOnly=TRUE))
xtable(cbind(myTable(marResultsLoad,0.2,c(1,2,3,5,6,8,10), ciOnly=TRUE), myTable(j2rResultsLoad,0.1,c(1,2,3,5,6,8,10), ciOnly=TRUE)),
       digits=c(3,3,2,3,2))