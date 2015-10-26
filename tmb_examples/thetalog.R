Y<-scan('thetalog.dat', skip=3, quiet=TRUE)

library(TMB)
compile("thetalog.cpp")
dyn.load(dynlib("thetalog"))
data <- list(Y=Y)
parameters <- list(
  X=data$Y*0,
  logr0=0,
  logtheta=0,
  logK=6,
  logQ=0,
  logR=0
  )
obj <- MakeADFun(data,parameters,random="X",DLL="thetalog")
newtonOption(obj,smartsearch=FALSE)
obj$fn()
obj$gr()
system.time(opt <- nlminb(obj$par,obj$fn,obj$gr))
rep <- sdreport(obj)
rep

