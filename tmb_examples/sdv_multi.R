require(TMB)

# Read data
source("sdv_multi_data.R")

library(TMB)
compile("sdv_multi.cpp")
dyn.load(dynlib("sdv_multi"))
obj <- MakeADFun(data=
                 list(n=n,p=p,y=t(y)),
                 parameters=list(
		   phi=rep(0.97,p),
		   log_sigma=rep(-1.7,p),
		   mu_x=rep(-0.5,p),
		   off_diag_x=rep(0.0,p),
    		   h=t(matrix(0.0,nrow=n,ncol=p))
                   ),
                 random=c("h")
#		   map=list(off_diag_x=factor(rep(NA,3)))
		 )
#obj$control <- list(trace=1,parscale=rep(1,12)*1e-2,REPORT=1,reltol=1e-12,maxit=100)

obj$hessian <- F
#ttt=obj$fn()
##Rprof();opt <- do.call("optim",obj);Rprof(NULL)
newtonOption(smartsearch=TRUE)
##print(system.time(opt <- do.call("optim",obj)))
ttt=system.time(opt<-nlminb(obj$par,obj$fn,obj$gr,
         lower=c(rep(-.99,3),rep(-3.0,3),rep(-3.0,3),rep(-5.0,3)),
         upper=c(rep(.99,3),rep(3.0,3),rep(3.0,3),rep(5.0,3))))
##obj$gr()
#c(phi1,phi2)
#f(opt$par)
rep <- sdreport(obj)
