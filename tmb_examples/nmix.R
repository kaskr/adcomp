source("tools/readdat.R")
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

# There is no direct equivalent of "phases" in TMB, but it can be achieved via the 
# "map" argument to MakeADFun. "phases" means keeping a subset of the parameters
# fixed when the model is first fit.

## Phase 1: u and log_sigma fixed
map <- list(u=factor(rep(NA,data$nG)),log_sigma=factor(NA))    # Say which parameters are fixed
obj <- MakeADFun(data,parameters,map=map,DLL="nmix")           # Note "map" argument
opt <- nlminb(obj$par,obj$fn,obj$gr)                           # Fit model
pl <- obj$env$parList(opt$par)                            # Parameter estimate after phase 1
## Phase 2: all parameters active
obj <- MakeADFun(data,pl,random="u",DLL="nmix")           # Note "pl" and missing "map"
system.time( opt <- nlminb(obj$par,obj$fn,obj$gr) )       # Fit model again
rep <- sdreport(obj)
rep
