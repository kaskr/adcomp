library(TMB)
compile("laplace.cpp")
dyn.load(dynlib("laplace"))

## Read data
source("spatial_data.R")

## Euclidian distance matrix
dd <- as.matrix(dist(Z))

## Data and parameter objects for TMB
data <- list(y=y, X=X, dd=dd, niter=5)
parameters <- list(
    b         = c(0, 0),
    a         = 1.428571,
    log_sigma = -0.6931472
)

## Fit model
obj <- MakeADFun(data, parameters, DLL="laplace")
system.time(opt<-nlminb(obj$par,obj$fn,obj$gr))
rep <- sdreport(obj)
rep
