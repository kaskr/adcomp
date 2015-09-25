library(TMB)
compile("laplace.cpp")
dyn.load(dynlib("laplace"))

## Read data
source("spatial_data.R")
dd <- sqrt(outer(Z[,1],Z[,1],"-")^2 + outer(Z[,2],Z[,2],"-")^2)
obj <- MakeADFun(data=list(y=y,X=X,dd=dd,niter=5),
                 parameters=list(
                     b=c(0,0),
                     a=1.428571,
                     log_sigma=-0.6931472
                     )
		 )
system.time(opt<-nlminb(obj$par,obj$fn,obj$gr))
rep <- sdreport(obj)
rep
