rm(obj);gc()
set.seed(123)
n <- 100
ng <- 10
f <- gl(ng,n/ng)
t <- rnorm(n,mean=2,sd=5)
a <- rnorm(ng,mean=3,sd=4)
b <- rnorm(ng,mean=8,sd=7)
x <- a[f]*t+b[f]+rnorm(n,mean=0,sd=1)
if(FALSE){
  library(lattice);xyplot(x~t,group=f,type="l")
}
require(TMB)
dyn.load(dynlib("randomregression"))

## Construct AD function object
data <- list(group=f,x=x,t=t)
parameters <- list(
                   a=rep(1,nlevels(f)),
                   b=rep(1,nlevels(f)),
                   mu=rep(1,2),
                   sigma=rep(1,2),
                   sigma0=1
                   )
obj <- MakeADFun(data=data,
                 parameters=parameters,
                 random=c("a","b"),
                 type = c("ADFun", "Fun"),
                 DLL="randomregression"
                 )

obj$fn()
obj$gr()

## Estimate
obj$control <- list(trace=1)
obj$hessian <- TRUE
newtonOption(smartsearch=FALSE,maxit=5)
##Rprof()
system.time(opt <- do.call("optim",obj))
##Rprof(NULL)
opt$par

## First four should be close to:
c(mean(a),mean(b),sd(a),sd(b))

