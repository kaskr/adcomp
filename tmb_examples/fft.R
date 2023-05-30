library(TMB)

## Compile and load the model
compile("fft.cpp")
dyn.load(dynlib("fft"))

## Data and parameters
set.seed(1)
data <- list(x=0:99, y=rnorm(100))
parameters <- list(rho=.1)

## Make a function object
obj <- MakeADFun(data, parameters, DLL="fft")

## Call function minimizer
fit <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr
