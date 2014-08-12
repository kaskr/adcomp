source("readdat.R")

d=readadmb("../admb_examples/socatt/socatt.dat")

d$ngroup = d$ngroup[-length(d$ngroup)]

d$y=as.factor(d$y)

d$X=matrix(d$X,ncol=d$p,byrow=TRUE)

group = as.factor(rep(1:d$M,each=4))

library(TMB)

compile("socatt.cpp")
dyn.load(dynlib("socatt"))

data <- list(y=d$y,S=d$S,X=d$X,group=group)
parameters <- list(
  u=rep(0,nlevels(group)),
  b=rep(0,d$p),
  logsigma=1,
  tmpk=rep(0,d$S-1)
  )
obj <- MakeADFun(data,parameters,random="u",DLL="socatt")
obj$fn()
obj$gr()
system.time(opt <- do.call("optim",obj))
pl <- obj$env$parList() ## <-- List of predicted random effects
rep <- sdreport(obj)
rep
