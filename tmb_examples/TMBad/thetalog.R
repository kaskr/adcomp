## Demonstrate automatic HMM filter (sequential reduction) on 'thetalog'

library(TMB)
compile("thetalog.cpp")
dyn.load(dynlib("thetalog"))
## Read data
Y <- scan("thetalog.dat", skip=3, quiet=TRUE)
data <- list(Y=Y, flag=0)

## Parameter initial guess
parameters <- list(
  X = data$Y*0,
  logr0 = 0,
  logtheta = 0,
  logK = 6,
  logQ = 0,
  logR = 0
)

## Specify configuration for sequential reduction
integrate <- list("X" = list("SR", "continuous", seq(0, 10, length=100)))

## Construct object with 'X' integrated by sequential reduction
obj <- MakeADFun(data, parameters, random="X", integrate=integrate, DLL="thetalog")

## Fit model
system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))

## Turn on epsilon method to get posterior mean
data <- list(Y=Y, flag=1)
parameters$eps <- parameters$X * 0
##parameters$scale <- 1
obj <- MakeADFun(data, parameters, random="X", integrate=integrate, DLL="thetalog")
obj$par[1:5] <- opt$par
obj$fn(obj$par)
g <- obj$gr(obj$par)

plot(tail(as.vector(g),length(Y)))
## load("../thetalog.expected.RData")
## points(tail(summary(.results$`TMB::sdreport`)[,1], length(Y)), col="red", pch=".")
