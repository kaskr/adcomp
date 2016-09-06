# Scaled up version of the Orange Tree example (5,000 latent random variables)

library(TMB)
compile("orange_big.cpp")
dyn.load(dynlib("orange_big"))

# Read data
source("orange_data.R")
Mmultiply <- data_orange$M * data_orange$multiply

parameters <- list(
    beta        = c(0,0,0),
    log_sigma   = 1,
    log_sigma_u = 2,
    u           = rep(0, Mmultiply)
)
obj <- MakeADFun(data_orange, parameters, random=c("u"), DLL="orange_big")

system.time(
    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  lower = c(-10.0, -10, -10, -5.0, -5.0),
                  upper = c( 10.0,  10,  10,  5.0,  5.0) )
)
rep <- sdreport(obj)
rep
