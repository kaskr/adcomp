## Simulate data
set.seed(123)
n <- 10000
mu <- 3
sigma <- 1.5
y <- rnorm(n, mu, sigma)

## data and parameters
data <- list(y=y)
parameters <- list(mu=0, logSigma = log(1), s=numeric(n))

## Compile C++ code and load into R
library(TMB)
compile("spa_gauss.cpp")
dyn.load(dynlib("spa_gauss"))

## create adfun, set s="random" for SPA inner problem
obj <- MakeADFun(data, parameters, random="s", DLL="spa_gauss",
                 intern=TRUE, inner.control=list(SPA=TRUE, sparse=TRUE))

## optimize
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)
## FIXME: waiting for delta method intern=TRUE
## rbind(summary(rep, "fixed", p.value = TRUE), summary(rep, "report", p.value = TRUE))

## Check values of nll vs nlspa
obj$fn(opt$par) 
-sum(dnorm(y, opt$par[1], exp(opt$par[2]), log=TRUE))

## Check values of dnll vs dnlspa
obj$gr(opt$par)
c(
    sum((opt$par[1]-y)/exp(opt$par[2])^2), # dnll / dmu
    -sum(exp(-2*opt$par[2])*(y-opt$par[1])^2) + n # dnll / dlogsigma
)

## check s = normal sp and parameters
##   arg min K - sx
##   K' - x = 0
##   K' = mu + sigma^2 s = x
##   s = (x - mu)/sigma^2
## FIXME: waiting for delta method intern=TRUE
## plot(summary(rep, "random")[,"Estimate"], (y-opt$par[1])/(exp(opt$par[2]))^2,
##      main="Numerical versus theoretical saddlepoints",
##      xlab="Numerical", ylab="Theoretical")
