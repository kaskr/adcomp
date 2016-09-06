library(TMB)
compile("tweedie.cpp")
dyn.load(dynlib("tweedie"))

## For reproducibility across machines:
if(TRUE) {
    data <- dget("tweedie_data.R")
} else {
    ## Simulate data
    set.seed(1001)
    data <- list( y = tweedie:::rtweedie(1000, mu=2, phi=2, p=1.5) )
}

## Parameter initial guess
parameters <- list(mu=1.1, phi=1.1, p=1.1)

## Fit model
model <- MakeADFun(data, parameters, DLL="tweedie")
model$fn()
model$gr()
system.time( opt <- nlminb(model$par, model$fn, model$gr) )
sdreport(model)
