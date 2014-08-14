# Scaled up version of the Orange Tree example (5,000 latent random variables)

require(TMB)

# Read data
source("data_orange.R")

compile("orange_big.cpp")
dyn.load(dynlib("orange_big"))
Mmultiply = data_orange$M*data_orange$multiply

obj <- MakeADFun(data=data_orange,
                 parameters=list(
		   beta=c(0,0,0),
		   log_sigma=1,
		   log_sigma_u=2,
                   u = rep(0,Mmultiply)),	
                 random=c("u")
		 )
obj$control <- list(trace=1,parscale=c(1,1,1,1,1)*1e-2,REPORT=1,reltol=1e-12,maxit=100)
obj$hessian <- F
newtonOption(smartsearch=TRUE)
tid = system.time(opt<-nlminb(obj$par,obj$fn,obj$gr,lower=c(-10.0,-10,-10,-5,-5),upper=c(10.0,10,10,5.0,5.0)))
rep <- sdreport(obj)
