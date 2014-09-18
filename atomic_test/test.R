library(TMB)
library(numDeriv)
##compile("test.cpp","-O0 -g -Woverloaded-virtual")
compile("test.cpp")
dyn.load(dynlib("test"))
##config(optimize.instantly=0)
x <- (2:9) / 10

## =============== Test pnorm
obj <- MakeADFun(data=list(a=0),parameters=list(x=x),DLL="test")
obj$fn(x)
sum(pnorm(x)) ## Check
obj$gr(x)
dnorm(x) ## Check

## 2nd order test
obj <- MakeADFun(data=list(a=0),parameters=list(x=x),DLL="test",random="x")
obj$env$spHess(x)
library(numDeriv)
diag(hessian(obj$env$f,x))

## =============== Test qnorm
obj <- MakeADFun(data=list(a=1),parameters=list(x=x),DLL="test")
obj$fn(x)
sum(qnorm(x)) ## Check
obj$gr(x)
1/(dnorm(qnorm(x))) ## Check

## 2nd order test
obj <- MakeADFun(data=list(a=1),parameters=list(x=x),DLL="test",random="x")
obj$env$spHess(x)
diag(hessian(obj$env$f,x))

## =============== Test incomplete gamma
incpl_gamma <- function(x)pgamma(x[1],x[2])*gamma(x[2])
obj <- MakeADFun(data=list(a=2),parameters=list(x=x),DLL="test")
obj$fn(x)
incpl_gamma(x[1:2]) ## Check
obj$gr(x)
grad(incpl_gamma,x[1:2]) ## Check

## 2nd order test
obj <- MakeADFun(data=list(a=2),parameters=list(x=x),DLL="test",random="x")
obj$env$spHess(x)
hessian(incpl_gamma,x[1:2]) ## Check

## =============== Test matrix multiply
f <- function(x){
    n <- length(x)/2
    m1 <- matrix(x[1:n],sqrt(n))
    m2 <- matrix(x[-(1:n)],sqrt(n))
    sum(m1%*%m2)
}
obj <- MakeADFun(data=list(a=3),parameters=list(x=x),DLL="test")
obj$fn(x)
f(x) ## Check
obj$gr(x)
grad(f,x) ## Check

## 2nd order test
obj <- MakeADFun(data=list(a=3),parameters=list(x=x),DLL="test",random="x")
obj$env$spHess(x)
round(hessian(f,x),8) ## Check

## =============== Test matrix inverse
f <- function(x){
    n <- length(x)
    m <- matrix(x,sqrt(n),sqrt(n))
    sum(solve(m))
}
set.seed(123);x <- rnorm(9)
obj <- MakeADFun(data=list(a=4),parameters=list(x=x),DLL="test")
obj$fn(x)
f(x) ## Check
obj$gr(x)
grad(f,x) ## Check

## 2nd order test
obj <- MakeADFun(data=list(a=4),parameters=list(x=x),DLL="test",random="x")
obj$env$spHess(x)
hessian(f,x) ## Check

## =============== Test matrix log-determinant
f <- function(x){
    n <- length(x)
    m <- matrix(x,sqrt(n),sqrt(n))
    determinant(m)$mod
}
set.seed(123);x <- matrix(rnorm(9),3);x <- x%*%t(x);x <- as.vector(x) ## PD only !
obj <- MakeADFun(data=list(a=5),parameters=list(x=x),DLL="test")
obj$fn(x)
f(x) ## Check
obj$gr(x)
grad(f,x) ## Check

## =============== Test inv incomplete gamma
x[1:2] <- c(.5,2)
inv_incpl_gamma <- function(x)qgamma(x[1]/gamma(x[2]),x[2])
obj <- MakeADFun(data=list(a=6),parameters=list(x=x),DLL="test")
obj$fn(x)
inv_incpl_gamma(x[1:2]) ## Check
obj$gr(x)
grad(inv_incpl_gamma,x[1:2]) ## Check

## =============== Test lgamma
obj <- MakeADFun(data=list(a=7),parameters=list(x=x),DLL="test")
obj$fn(x)
lgamma(x[1]) ## Check
obj$gr(x)
grad(lgamma,x[1]) ## Check

## =============== Test pgamma
x[1:3] <- 1:3
f <- function(x)pgamma(x[1],x[2],scale=x[3])
obj <- MakeADFun(data=list(a=8),parameters=list(x=x),DLL="test")
obj$fn(x)
f(x[1:3]) ## Check
obj$gr(x)
grad(f,x[1:3]) ## Check

## =============== Test qgamma
x[1:3] <- c(.7,2,3)
f <- function(x)pgamma(x[1],x[2],scale=x[3])
obj <- MakeADFun(data=list(a=8),parameters=list(x=x),DLL="test")
obj$fn(x)
f(x[1:3]) ## Check
obj$gr(x)
grad(f,x[1:3]) ## Check
