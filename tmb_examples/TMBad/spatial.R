library(TMB)
compile("spatial.cpp")
dyn.load(dynlib("spatial"))

## Read data
source("spatial_data.R")

## Euclidian distance matrix
dd <- as.matrix(dist(Z))

obj <- MakeADFun(data = list(n=100, y=y, X=X, dd=dd),
                 parameters=list(
                     b = c(0,0),
                     a = 1.428571,
                     log_sigma = -0.6931472,
                     u = rep(0,n)),
                 DLL = "spatial",
                 random = "u" , intern=TRUE,
                 inner.control=list(sparse=FALSE)
		 )

system.time(
    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  lower=c(-100.0, -100.0, 0.01, -3.0),
                  upper=c( 100.0,  100.0, 3.00,  3.0) )
)

rep <- sdreport(obj)
rep
