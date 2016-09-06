library(TMB)
compile("sdv_multi.cpp")
dyn.load(dynlib("sdv_multi"))

# Read data
source("sdv_multi_data.R")
data <- list(n=n, p=p, y=t(y))

## Parameter initial guess
parameters <- list(
    phi        = rep(0.97,p),
    log_sigma  = rep(-1.7,p),
    mu_x       = rep(-0.5,p),
    off_diag_x = rep(0.0,p),
    h          = t(matrix(0.0,nrow=n,ncol=p))
)

## Fit model
obj <- MakeADFun(data, parameters, random="h", DLL="sdv_multi")
system.time(
    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  lower = c(rep(-.99,3), rep(-3.0,3), rep(-3.0,3), rep(-5.0,3)),
                  upper = c(rep( .99,3), rep( 3.0,3), rep( 3.0,3), rep( 5.0,3)) )
)
rep <- sdreport(obj)
rep
