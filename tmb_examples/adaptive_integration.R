require(TMB)
compile("adaptive_integration.cpp")
dyn.load(dynlib("adaptive_integration"))

## Simulate data
set.seed(123)
ndat <- 1000
c1 <- rnorm(ndat, sd=1)
c2 <- rnorm(ndat, sd=1)
ngroup <- 11
c3 <- cut(runif(ndat), ngroup)
eps <- rnorm(ndat, sd=.5)
p <- plogis(c1 + c2 + seq(-1,1,length=ngroup)[c3] + eps)
n <- rep(10, ndat)
x <- rbinom(ndat, size=n, prob=p)
A <- model.matrix(~ c1 + c2 + c3 - 1)

## Data and parameters for TMB
data <- list(x=x, n=n, A=A)
parameters <- list(b = rep(0,ncol(A)), logsd=0)

## Fit model
model <- MakeADFun(data, parameters, DLL="adaptive_integration")
system.time( fit <- nlminb(model$par, model$fn, model$gr) )
rep <- sdreport(model)
print(rep)
