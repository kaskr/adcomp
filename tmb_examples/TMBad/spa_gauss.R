## Simulate data
set.seed(123)
n <- 1000
mu <- 3
sigma <- 1.5
y <- rnorm(n, mu, sigma)

## data and parameters
data <- list(y=y,s=numeric(n))
parameters <- list(mu=0, logSigma = log(1))

## Compile C++ code and load into R
library(TMB)
compile("spa_gauss.cpp", framework="TMBad")
dyn.load(dynlib("spa_gauss"))

## create adfun (SPA inner problem option set in C++ template)
obj <- MakeADFun(data, parameters, DLL="spa_gauss")

## optimize
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)
summary(rep, p.value=TRUE)

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
s_hat <- obj$report()$s
plot(s_hat, (y-opt$par[1])/(exp(opt$par[2]))^2,
     main="Numerical versus theoretical saddlepoints",
     xlab="Numerical", ylab="Theoretical")
