library(TMB)

## Compile and load the model
compile("solver.cpp")
dyn.load(dynlib("solver"))

## Data and parameters
n <- 5
m <- diag(n) + 5

x <- 1:n
data <- list(m=m)
parameters <- list(x=x)

## Make a function object
obj <- MakeADFun(data, parameters, DLL="solver")

## This is what tape looks like:
TMB:::tape_print(obj$env$ADFun)

## We can differentiate solution to any order!
obj$fn()
obj$gr()
obj$he()

## Compare with known solution
library(numDeriv)
f <- function(x){sum(x*x)}
f(x)
grad(f,x)
round(hessian(f,x),10)
