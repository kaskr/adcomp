## Generate AR1 variables
set.seed(123)
n <- 1000
phi <- .6
u <- numeric(n)
u[0] = rnorm(1)
for(i in 2:n){
    u[i] = phi * u[i-1] + rnorm(1, sd = sqrt(1 - phi^2))
}

## Transform to gamma
shape <- 2
scale <- 3
x <- qgamma(pnorm(u),shape=shape,scale=scale)

## Samples with noise
sd <- 2
y <- x + rnorm(n,sd=sd)

library(TMB)
compile("transform.cpp")
dyn.load(dynlib("transform"))

data=list(y=y)
parameters=list(phi=0,shape=1,scale=1,sd=1,u=u*0)
obj <- MakeADFun(
    data=data,parameters=parameters,
    DLL="transform",
    random="u"
    )
obj$fn()
opt <- nlminb(obj$par,obj$fn,obj$gr)
sdr <- sdreport(obj)
summary(sdr,"fixed")

