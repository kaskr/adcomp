source("sam_data.R")

library(TMB)
compile("sam.cpp")
dyn.load(dynlib("sam"))

parameters <- list(
  logFpar      = c(-11, -10, -10, -9, -9),
  logQpow      = numeric(0),
  logSdLogFsta = rep(-0.693147,max(data$keyVarF)+1),
  logSdLogN    = c(0.356675, -0.356675),
  logSdLogObs  = c(-0.356675, -0.356675, -0.356675, -0.356675, -0.356675),
  rec_loga     = 1,
  rec_logb     = -12,
  rho          = 0.5,
  logScale     = numeric(data$noScaledYears),
  logScaleSSB  = if(any(data$fleetTypes %in% c(3,4))) {numeric(1)} else {numeric(0)},
  logPowSSB    = if(any(data$fleetTypes == 4))        {numeric(1)} else {numeric(0)},
  logSdSSB     = if(any(data$fleetTypes %in% c(3,4))) {numeric(1)} else {numeric(0)},
  U = matrix(0, nrow=max(data$keyLogFsta)+1 + data$maxAge-data$minAge+1, ncol=data$noYears)
)
data$nlogF = max(data$keyLogFsta)+1
data$nlogN = data$maxAge-data$minAge+1

obj <- MakeADFun(data, parameters, random=c("U"), DLL="sam", intern=TRUE)
lower <- obj$par*0-Inf
upper <- obj$par*0+Inf
lower["rho"] <- 0.01
upper["rho"] <- 0.99

system.time(opt <- nlminb(obj$par, obj$fn, obj$gr, lower=lower, upper=upper))
hessian <- obj$he(opt$par)
rep <- sdreport(obj, hessian=hessian)
rep
