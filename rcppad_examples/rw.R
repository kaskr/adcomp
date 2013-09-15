set.seed(1);
library(MASS)

simdata <- function(){
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
stateDim=3
timeSteps=100
simdata()
## make ADMB data file
## cat(timeSteps,"\n",file="mvrw.dat");
## cat(stateDim,"\n",file="mvrw.dat",append=TRUE);
## write.table(obs,file="mvrw.dat",append=TRUE,row.names=FALSE,col.names=FALSE)
library(RcppAD)
compile("rw.cpp")
compile("rw_parallel.cpp")
dyn.load("rw.so")
dyn.load("rw_parallel.so")
data <- list(obs=t(obs))
parameters <- list(
  u=data$obs*0,
  transf_rho=0.1,
  logsds=sds*0,
  logsdObs=sdObs*0
  )
newtonOption(smartsearch=FALSE)
obj <- MakeADFun(data,parameters,random="^u",DLL="rw")
obj$hessian <- TRUE

obj$fn()
obj$gr()
system.time(opt <- do.call("optim",obj))
pl <- obj$env$parList() ## <-- List of predicted random effects
matpoints(t(pl$u),type="l",col="blue",lwd=2)
summary(as.mle(opt))


## Permormance analysis in case 1000 timesteps and 10 states
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
obj <- MakeADFun(data,parameters,random="^u",DLL="rw_parallel")
ben <- benchmarkParallel(obj,n=1000,cores=1:8)
ben
plot(ben)
