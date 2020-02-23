library(TMB)

## Compile and load the model
compile("solver.cpp")
dyn.load(dynlib("solver"))

## Data and parameters
n <- 5
m <- diag(n) + 5

x <- rep(1, n)
data <- list(m=m, trace=1)
parameters <- list(x=x)

## Make a function object
obj <- MakeADFun(data, parameters, DLL="solver")

## This is what tape looks like:
TMB:::tape_print(obj$env$ADFun)

## We can differentiate solution to any order!
obj$fn()
obj$gr()
obj$he()

## Check derivatives
library(numDeriv)
p <- obj$par
g <- grad(obj$fn, p)
h <- jacobian(obj$gr, p)
range(obj$gr(p)-g)
range(obj$he(p)-h)

nlminb(p, obj$fn, obj$gr)
