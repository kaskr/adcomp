source("readdat.R")
data <- readadmb()
data$IDfac <- local(factor(rep(1:length(nID),nID)),data)
length(data$IDfac)==length(data$ID)
data$IDind <- matrix(data$IDind,data$R,data$T,byrow=TRUE)
data$y <- matrix(data$y,data$R,data$T,byrow=TRUE)
data$IDind <- data$IDind-1
data$ID <- factor(data$ID)
library(TMB)
compile("nmix.cpp")
parameters <- list(
                   log_lambda=4,
                   p0=0,
                   p1=0,
                   log_sigma=0,
                   u=rep(0,data$nG)
)
dyn.load(dynlib("nmix"))
## Phase 1
map <- list(u=factor(rep(NA,data$nG)),log_sigma=factor(NA))
obj <- MakeADFun(data,parameters,map=map,DLL="nmix")
opt <- nlminb(obj$par,obj$fn,obj$gr)
pl <- obj$env$parList(opt$par) ## Parameter estimate after phase 1
## Phase 2
obj <- MakeADFun(data,pl,random="u",DLL="nmix")
system.time( opt <- nlminb(obj$par,obj$fn,obj$gr) )
rep <- sdreport(obj)
rep
