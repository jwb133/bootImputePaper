#this file analyses the simulated datasets in the first simulation study in the paper

#the code below is for running on a high performance cluster using the SLURM scheduling system
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
# coerce the value to an integer
batch <- as.numeric(slurm_arrayid)

#find seed for this run
library(parallel)
RNGkind("L'Ecuyer-CMRG")
set.seed(69012365) # something
s <- .Random.seed
for (i in 1:batch) {
  s <- nextRNGStream(s)
}
.GlobalEnv$.Random.seed <- s

print(batch)

# read in simulated data
load(file=paste("./datasets/inputData_",batch,".RData",sep=""))
n <- dim(originalDataList[[1]])[1]

library(MASS)
library(mlmi)
library(mitools)


#function that imputes Y2 using J2R with or without posterior draw
j2rImp <- function(inputData,postDraw) {
  z <- inputData[,1]
  y <- inputData[,c(2,3)]
  #we first perform standard MAR imputation for control group
  controly <- data.matrix(y[z==0,])
  ymod <- lm(y.2~y.1, data=data.frame(controly))
  if (postDraw==FALSE) {
    outcomeModBeta <- coef(ymod)
    outcomeModResVar <- summary(ymod)$sigma^2
  } else {
    #posterior draw first
    beta <- ymod$coef
    sigmasq <- summary(ymod)$sigma^2
    varcov <- vcov(ymod)
    outcomeModResVar <- (sigmasq*ymod$df) / rchisq(1,ymod$df)
    covariance <- (outcomeModResVar/sigmasq)*vcov(ymod)
    outcomeModBeta = beta + MASS::mvrnorm(1, mu=rep(0,ncol(covariance)), Sigma=covariance)
  }
  impfitted <- outcomeModBeta[1] + outcomeModBeta[2]*controly[,1]
  controly[is.na(controly[,2])] <- impfitted[is.na(controly[,2])] + rnorm(sum(is.na(controly[,2])), mean=0, sd=outcomeModResVar^0.5)
  
  #now j2r for active patients
  #this uses regression coefficient from controls, marginal mean for y2 from control, and marginal mean for y1 from active
  activey <- data.matrix(y[z==1,])
  #control and active means at time 1
  if (postDraw==FALSE) {
    mu1Control <- mean(controly[,1])
    mu1Active <- mean(activey[,1])
  } else {
    SigmaSq1Control <- (dim(controly)[1]-1)*var(controly[,1])/rchisq(1,dim(controly)[1]-1)
    mu1Control <- rnorm(1,mean=mean(controly[,1]), sd=(SigmaSq1Control/dim(controly)[1])^0.5)
    SigmaSq1Active <- (dim(activey)[1]-1)*var(activey[,1])/rchisq(1,dim(activey)[1]-1)
    mu1Active <- rnorm(1,mean=mean(activey[,1]), sd=(SigmaSq1Active/dim(activey)[1])^0.5)
  }
  #control mean at time 2
  mu2Control <- outcomeModBeta[1] + outcomeModBeta[2]*mu1Control
  #fitted values for active patients under j2r
  impfitted <- mu2Control + outcomeModBeta[2]*(activey[,1]-mu1Active)
  activey[is.na(activey[,2])] <- impfitted[is.na(activey[,2])] + rnorm(sum(is.na(activey[,2])), mean=0, sd=outcomeModResVar^0.5)
  
  data.frame(rbind(cbind(z=0, controly), cbind(z=1, activey)))
}

bootMIvonHippel <- function(estsArray, B, M) {
  #fit one way model
  SSW <- sum((estsArray-rowMeans(estsArray))^2)
  SSB <- M*sum((rowMeans(estsArray)-mean(estsArray))^2)
  MSW <- SSW/(B*(M-1))
  MSB <- SSB/(B-1)
  resVar <- MSW
  randIntVar <- (MSB-MSW)/M
  if (randIntVar<0) {
    randIntVar <- 0
    resVar <- (SSW+SSB)/(B*M-1)
  }
  pointEstimate <- mean(estsArray)
  varEstimate <- (1+1/B)*randIntVar + resVar/(B*M)
  df <- (varEstimate^2)/((((B+1)/(B*M))^2*MSB^2 / (B-1)) + MSW^2/(B*M^2*(M-1)))
  c(pointEstimate,varEstimate, pointEstimate-qt(0.975,df=df)*varEstimate^0.5, pointEstimate+qt(0.975,df=df)*varEstimate^0.5)
}


imputeAnalyse <- function(mar) {
  
  #standard Rubin, MI boot pooled, MI boot Rubin
  for (i in 1:bsM) {
    if (mar==TRUE) {
      yimp <- normUniImp(originalData, y.2~z+y.1, M=1, pd=TRUE)
    } else {
      yimp <- j2rImp(originalData,postDraw=TRUE)
    }
    
    impDataMod <- impDataMod <- lm(y.2~z+y.1, data=yimp)
    impEsts[i] <- coef(impDataMod)[2]
    impVars[i] <- vcov(impDataMod)[2,2]
    
    for (j in 1:bsSamples) {
      #take bootstrap sample of yimp
      bsData <- yimp[sample(n, replace=TRUE),]
      impDataMod <- impDataMod <- lm(y.2~z+y.1, data=bsData)
      bsImpEsts[j,i] <- coef(impDataMod)[2]
    }
    impVars2[i] <- var(bsImpEsts[,i])
  }
  rubinResult <- summary(MIcombine(results=as.list(impEsts), variances=as.list(impVars), df.complete=n-3))
  rubin <- c(rubinResult$results, rubinResult$se^2, rubinResult$`(lower`, rubinResult$`upper)`)
  
  rubinResult <- summary(MIcombine(results=as.list(impEsts), variances=as.list(impVars2), df.complete=n-3))
  miBootRubin <- c(rubinResult$results, rubinResult$se^2, rubinResult$`(lower`, rubinResult$`upper)`)
  
  miBootPooledPercentile <- c(mean(impEsts), var(c(bsImpEsts)), quantile(c(bsImpEsts),0.025), quantile(c(bsImpEsts),0.975))
  miBootPooledNormal <- c(mean(impEsts), var(c(bsImpEsts)), mean(impEsts)-1.96*var(c(bsImpEsts))^0.5, 
                         mean(impEsts)+1.96*var(c(bsImpEsts))^0.5)
  
  ##################################################Boot MI methods####################################################
  #boot MI and boot MI pooled
  for (i in 1:bsSamples) {
    bsData <- originalData[sample(n, replace=TRUE),]
    
    for (j in 1:bsM) {
      if (mar==TRUE) {
        yimp <- normUniImp(bsData, y.2~z+y.1, M=1, pd=TRUE)
      } else {
        yimp <- j2rImp(bsData,postDraw=TRUE)
      }
      impDataMod <- impDataMod <- lm(y.2~z+y.1, data=yimp)
      bsImpEsts[i,j] <- coef(impDataMod)[2]
    }
  }
  bootMIpercentile <- c(rubin[1], var(rowMeans(bsImpEsts)), quantile(rowMeans(bsImpEsts),0.025),quantile(rowMeans(bsImpEsts),0.975))
  bootMI1percentile <- c(impEsts[1], var(bsImpEsts[,1]), quantile(bsImpEsts[,1],0.025),quantile(bsImpEsts[,1],0.975))
  
  bootMInormal <- c(rubin[1], var(rowMeans(bsImpEsts)), rubin[1]-1.96*var(rowMeans(bsImpEsts))^0.5,
                    rubin[1]+1.96*var(rowMeans(bsImpEsts))^0.5)
  bootMI1normal <- c(impEsts[1], var(bsImpEsts[,1]), impEsts[1]-1.96*var(bsImpEsts[,1])^0.5, impEsts[1]+1.96*var(bsImpEsts[,1])^0.5)
  
  bootMIpooledPercentile <- c(mean(bsImpEsts), var(bsImpEsts), quantile(c(bsImpEsts),0.025),quantile(c(bsImpEsts),0.975))
  bootMIpooledNormal <- c(mean(bsImpEsts), var(bsImpEsts), mean(bsImpEsts)-1.96*var(bsImpEsts)^0.5, mean(bsImpEsts)+1.96*var(bsImpEsts)^0.5)
  
  #von Hippel, first using all imputations within each bootstrap
  vonHippelAll <- bootMIvonHippel(bsImpEsts, bsSamples, bsM)
  
  #von Hippel only using M=2
  vonHippelM2 <- bootMIvonHippel(bsImpEsts[,c(1,2)], bsSamples, 2)
  
  rbind(rubin, miBootRubin, miBootPooledPercentile, miBootPooledNormal, bootMIpercentile, bootMI1percentile,
        bootMInormal, bootMI1normal, vonHippelAll, vonHippelM2=vonHippelM2)
}

#number of simulations per instance
nSim <- length(originalDataList)

#number of bootstraps
bsSamples <- 200
#number of imputations per bootstrap
bsM <- 10

bsImpEsts <- array(0, dim=c(bsSamples,bsM))

impEsts <- array(0, dim=bsM)
impVars <- array(0, dim=bsM)
impVars2 <- array(0, dim=bsM)

numMethods <- 10
marResults <- array(0, dim=c(nSim,numMethods,4))
j2rResults <- array(0, dim=c(nSim,numMethods,4))

for (sim in 1:nSim) {
  #print(sim)
  originalData <- originalDataList[[sim]]
  
  #first impute under MAR (congenial)
  marResults[sim,,] <- imputeAnalyse(mar=TRUE)
  
  #second impute using J2R (uncongenial)
  j2rResults[sim,,] <- imputeAnalyse(mar=FALSE)

}

save(marResults, j2rResults, file=(paste("./results/simRes_", batch, ".RData", sep="")))

