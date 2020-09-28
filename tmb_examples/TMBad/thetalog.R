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

## Fit model Laplace
if (FALSE) {
    obj <- MakeADFun(data, parameters, random="X", DLL="thetalog")
    newtonOption(obj, smartsearch=FALSE)
    obj$fn()
    obj$gr()
    system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))
    range(obj$env$last.par.best[1:length(Y)])
}

obj <- MakeADFun(data, parameters, DLL="thetalog")
random <- 1:length(Y)

#######
DLL <- "thetalog"
x <- seq(0,10,length=100)
w <- rep(diff(x)[1],length(x))
grid <- list(x=x,w=w)
## Construct HMM filter
con <- list(random_order = random, 
            method = "marginal_sr",
            grid = grid,
            mustWork = 1L,
            max_period_size = 1024L)
.Call("TransformADFunObject",
      obj$env$ADFun$ptr,
      con,
      PACKAGE = DLL)

## Remove random parameters from funciton objects
.Call("TransformADFunObject",
      obj$env$ADFun$ptr,
      list(method = "remove_random_parameters", 
           random_order = random,
           mustWork = 1L,
           max_period_size = 1024L), 
      PACKAGE = DLL)
attr(obj$env$ADFun$ptr, "par") <- attr(obj$env$ADFun$ptr, "par")[-random]
obj$env$par <- obj$env$par[-random]
obj$par <- obj$env$par
obj$env$random <- NULL

## Fit model
system.time(qw <- nlminb(obj$par, obj$fn, obj$gr)) ## 67.620
qw$evaluations

## Turn on epsilon method to get posterior mean
data <- list(Y=Y, flag=1)
parameters$eps <- parameters$X * 0
##parameters$scale <- 1
obj <- MakeADFun(data, parameters, DLL="thetalog")
random <- 1:length(Y)
.Call("TransformADFunObject",
      obj$env$ADFun$ptr,
      con,
      PACKAGE = DLL)
obj$par[-(c(1:length(Y)))][1:5] <- qw$par
obj$fn(obj$par)
g <- obj$gr(obj$par)

plot(tail(as.vector(g),length(Y)))
## load("../thetalog.expected.RData")
## points(tail(summary(.results$`TMB::sdreport`)[,1], length(Y)), col="red", pch=".")
