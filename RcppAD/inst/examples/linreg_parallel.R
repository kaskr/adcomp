library(RcppAD)
x <- seq(0,10,length=50001)
data <- list(Y=rnorm(length(x))+x,x=x)
parameters <- list(a=0,b=0,logSigma=0)
dyn.load("linreg_parallel.so")
obj <- MakeADFun(data,parameters,DLL="linreg_parallel")
obj$hessian <- TRUE
opt <- do.call("optim",obj)
summary(as.mle(opt))
