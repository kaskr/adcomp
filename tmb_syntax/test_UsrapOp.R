library(TMB)

## Compile and load the model
compile("test_UsrapOp.cpp")
dyn.load(dynlib("test_UsrapOp"))

## Data and parameters
n <- 400
m <- 50
set.seed(1)
data <- list(A=matrix(rnorm(n*m), n, m), y=rnorm(n))
parameters <- list(u=numeric(m), sd=1)

## Make a function object
obj <- MakeADFun(data, parameters, DLL="test_UsrapOp", random="u")

obj$fn()
obj$gr() ## Crash
