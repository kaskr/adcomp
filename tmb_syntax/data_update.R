## Demonstrate how to change data without re-taping
library(TMB)

## Compile and load example
compile("data_update.cpp")
dyn.load(dynlib("data_update"))

## Generate random objects
set.seed(123)
a <- rep(0, 2)
b <- array(0, c(2,3,4) )
c <- matrix(0,2,2)
d <- 1
x <- rnorm(1)
data <- list(a=a,b=b,c=c,d=d)
parameters <- list(x=x)

## Check that data objects can be changed from R
obj <- MakeADFun(data=data,parameters=parameters,DLL="data_update")

## Function to check value
check <- function(x) {
    res.true <- -dnorm(x, 0, sum(unlist(obj$env$data)), TRUE)
    res.obs  <- obj$fn(x)
    abs.diff <- abs(res.true - res.obs)
    stopifnot( abs.diff < 1e-10 )
    data.frame(res.true, res.obs, abs.diff)
}

## No change made yet:
check(x=1)
check(x=2)

## Change 'a' without re-taping
obj$env$data$a <- obj$env$data$a + .1
check(x=1)
check(x=2)

## Change 'b' without re-taping
obj$env$data$b <- obj$env$data$b + .01
check(x=1)
check(x=2)

## Change 'c' without re-taping
obj$env$data$c <- obj$env$data$c + .05
check(x=1)
check(x=2)

## Change 'd' without re-taping
obj$env$data$d <- obj$env$data$d + .1
check(x=1)
check(x=2)
