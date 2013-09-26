require(RcppAD)

#Modified to read files that I can run with ADMB

# Read data
source("spatial_data.R")

dd = sqrt(outer(Z[,1],Z[,1],"-")^2 + outer(Z[,2],Z[,2],"-")^2)	

library(RcppAD)
compile("spatial.cpp")
dyn.load("spatial.so")
obj <- MakeADFun(data=list(n=100,y=y,X=X,dd=dd),
                 parameters=list(
		   b=c(0,0),
                   a=1.428571,
		   log_sigma=-0.6931472,
                   u = rep(0,n)),	
                 random=c("u")
		 )
obj$control <- list(trace=1,parscale=c(1,1,1,1)*1e-2,REPORT=1,reltol=1e-12,maxit=100)

obj$hessian <- F
##Rprof();opt <- do.call("optim",obj);Rprof(NULL)
newtonOption(smartsearch=TRUE)
##print(system.time(opt <- do.call("optim",obj)))
system.time(opt<-nlminb(obj$par,obj$fn,obj$gr,lower=c(-100.0,-100.0,0.01,-3.0),upper=c(100,100,3.0,3.0)))
##obj$gr()
#c(phi1,phi2)
#f(opt$par)
