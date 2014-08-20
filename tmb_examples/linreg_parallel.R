library(TMB)
set.seed(123)
x <- seq(0,10,length=50001)
data <- list(Y=rnorm(length(x))+x,x=x)
parameters <- list(a=0,b=0,logSigma=0)
dyn.load(dynlib("linreg_parallel"))
obj <- MakeADFun(data,parameters,DLL="linreg_parallel")
obj$hessian <- TRUE
opt <- do.call("optim",obj)
