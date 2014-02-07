Y<-scan('thetalog.dat', skip=3, quiet=TRUE)

library(TMB)
compile("thetalog.cpp")
dyn.load("thetalog.so")
data <- list(Y=Y)
parameters <- list(
  X=data$Y*0,
  logr0=0,
  logtheta=0,
  logK=6,
  logQ=0,
  logR=0
  )
newtonOption(smartsearch=FALSE)
obj <- MakeADFun(data,parameters,random="^X",DLL="thetalog")
obj$hessian <- TRUE
#obj$control$reltol<-1e-12

obj$fn()
obj$gr()
system.time(opt <- nlminb(obj$par,obj$fn,obj$gr))

rep <- sdreport(obj)
rep

