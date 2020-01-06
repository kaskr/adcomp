library(TMB)
compile("compois.cpp")
dyn.load(dynlib("compois"))

set.seed(123)
nu <- .1
mode <- 10
domain <- 0:100
prob <- dpois(domain, lambda=mode)^nu; prob <- prob / sum(prob)
sum(prob * domain) ## mean
x <- sample(domain, size=1e4, replace=TRUE, prob = prob)

## TMB data
data <- list( x = x )
parameters <- list( logmu = 0, lognu = 0 )

## Parameterization through the mode
data$parameterization <- "mode"
obj <- MakeADFun(data, parameters, DLL="compois")
system.time(fit.mode <- nlminb(obj$par, obj$fn, obj$gr, obj$he))
rep.mode <- sdreport(obj)

## Parameterization through the mean
data$parameterization <- "mean"
obj <- MakeADFun(data, parameters, DLL="compois")
system.time(fit.mean <- nlminb(obj$par, obj$fn, obj$gr, obj$he))
rep.mean <- sdreport(obj)

summary(rep.mode, "report")
summary(rep.mean, "report")
