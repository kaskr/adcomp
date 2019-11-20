## Check correctness of 'autdiff' namespace
library(TMB)

## Compile and load the model
compile("check_autodiff.cpp")
dyn.load(dynlib("check_autodiff"))

## Data and parameters
data <- list()
parameters <- list(theta=1:2)

## Make a function object
data$select <- 1
obj <- MakeADFun(data, parameters, DLL="check_autodiff", ADreport=TRUE, silent=TRUE)

## Check autodiff::gradient and autodiff::hessian
stopifnot( all( obj$report(obj$par)$g == obj$gr(obj$par) ) )
stopifnot( identical( obj$report(obj$par)$h, obj$he(obj$par) ) )

## Make a function object
data$select <- 2
obj <- MakeADFun(data, parameters, DLL="check_autodiff", ADreport=TRUE, silent=TRUE)
stopifnot( all( obj$gr(obj$par) == diag(exp(obj$par)) ) )

## Make new function object
data$select <- 3
obj <- MakeADFun(data, parameters, DLL="check_autodiff", ADreport=TRUE, silent=TRUE)
stopifnot( all( obj$report(obj$par)$j == obj$gr(obj$par) ) )
