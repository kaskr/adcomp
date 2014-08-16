library(TMB)
#compile("longlinreg.cpp")
dyn.load(dynlib("longlinreg"))
set.seed(123)
nobs<-1000000
x<-seq(0,10, length=nobs)
data <- list(Y=2*x+1+rnorm(nobs),x=x)
parameters <- list(a=0,b=0,logSigma=0)
obj <- MakeADFun(data,parameters,DLL="longlinreg")
system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))
opt
obj$he()    ## <-- Analytical hessian
sdreport(obj)
