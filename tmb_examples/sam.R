source("sam_data.R")

library(TMB)
compile("sam.cpp")
dyn.load("sam.so")

parameters <- list(
  logFpar=c(-11, -10, -10, -9, -9),
  logQpow=numeric(0),
  logSdLogFsta=rep(-0.693147,max(data$keyVarF)+1),
  logSdLogN=c(0.356675, -0.356675),
  logSdLogObs=c(-0.356675, -0.356675, -0.356675, -0.356675, -0.356675),
  rec_loga=1,
  rec_logb=-12,
  logit_rho=0.55, 
  logScale=numeric(data$noScaledYears),  
  logScaleSSB=if(any(data$fleetTypes%in%c(3,4))){numeric(1)}else{numeric(0)},
  logPowSSB=if(any(data$fleetTypes==4)){numeric(1)}else{numeric(0)},
  logSdSSB=if(any(data$fleetTypes%in%c(3,4))){numeric(1)}else{numeric(0)},
  logF=matrix(0, nrow=max(data$keyLogFsta)+1,ncol=data$noYears),
  logN=matrix(0, nrow=data$maxAge-data$minAge+1, ncol=data$noYears)
)

obj <- MakeADFun(data,parameters,random=c("logN","logF"),DLL="sam")
system.time(opt<-nlminb(obj$par,obj$fn,obj$gr))
rep<-sdreport(obj)
rep
