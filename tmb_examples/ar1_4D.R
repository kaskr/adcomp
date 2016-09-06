require(TMB)
compile("ar1_4D.cpp")
dyn.load(dynlib("ar1_4D"))

set.seed(123)
n <- 8 ## Size of problem = n*n

## ======================= Simulate separable 2D GMRF 
## - With exponential correlation in both directions
## - phi1 = 1-lag correlation in 1st direction
## - phi2 = 1-lag correlation in 2nd direction
ar1corr <- function(n,phi){
  phi^abs(outer(1:n,1:n,"-"))
}
simgmrf4 <- function(n,phi){
  dim <- c(n,n,n,n)
  u <- array(rnorm(prod(dim)),dim)
  L <- t(chol(ar1corr(n,phi)))
  ar2mat <- function(x)matrix(x,nrow(x))
  for(i in 1:4){
    u[] <- L%*%ar2mat(u)
    u <- aperm(u,c(4,1,2,3))
  }
  u
}

## ======================= Simulate data
phi=exp(-1/(.2*n)) ## Correlation range=20% of grid size second dimension
eta <- simgmrf4(n,phi)
N <- rpois(length(eta),exp(eta))

## ======================= Parameterization of phi
f <- function(x) 2/(1 + exp(-2 * x)) - 1
invf <- function(y) -0.5 * log(2/(y + 1) - 1)

## ======================= Fit model
obj <- MakeADFun(data=list(N=N),
                 parameters=list(
                   eta=array(0,c(n,n,n,n)),
                   transf_phi=invf(0.5)
                   ),
                 random=c("eta"),
                 DLL="ar1_4D"
                 )
runSymbolicAnalysis(obj)
obj$control <- list(trace=1,parscale=c(1)*1e-2,REPORT=1,reltol=1e-12)
newtonOption(obj, smartsearch=FALSE)
system.time(opt <- do.call("optim",obj))
rep <- sdreport(obj)
rep
