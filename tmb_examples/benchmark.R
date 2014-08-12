library(MASS)

simdata <- function(){
  set.seed(1);
  local({
    rho=0.9
    sds=seq(0.5,2,length=stateDim)
    sdObs=rep(1,stateDim);
    corrMat=matrix(0.0,stateDim,stateDim)
    for(i in 1:stateDim){
      for(j in 1:stateDim){
        corrMat[i,j] = rho^abs(i-j)
      }
    }
    Sigma=corrMat*(sds %o% sds)
    d=matrix(NA,timeSteps,stateDim)
    obs=d;
    ##init state
    d[1,] = rnorm(stateDim);
    i=1;
    obs[i,] = d[i,] + rnorm(stateDim,rep(0,stateDim),sdObs)
    for(i in 2:timeSteps){
      d[i,] = d[i-1,] + mvrnorm(1,rep(0,stateDim),Sigma=Sigma)
      obs[i,] = d[i,] + rnorm(stateDim,rep(0,stateDim),sdObs)
    }
    matplot(d,type="l")
    matpoints(obs);
  },.GlobalEnv)
}

## Performance analysis in case 1000 timesteps and 10 states
library(TMB)
stateDim=10
timeSteps=1000
simdata()
data <- list(obs=t(obs))
parameters <- list(
  u=data$obs*0,
  transf_rho=0.1,
  logsds=sds*0,
  logsdObs=sdObs*0
  )
newtonOption(smartsearch=FALSE)
##compile("rw_parallel.cpp")
dyn.load(dynlib("rw_parallel"))
obj <- MakeADFun(data,parameters,random="u",DLL="rw_parallel")
ben <- benchmarkParallel(obj,n=1000,cores=1:8)
ben
pdf("../slides/results/scalability.pdf")
plot(ben,type="b",ylim=c(0,1.1),las=1)
plot(function(x)1/x,1,10,col="grey",add=TRUE)
legend("bottomleft",legend=rownames(ben),col=1:4,lty=1:4)
title("Scalability: stateDim=10 timeSteps=1000")
dev.off()
