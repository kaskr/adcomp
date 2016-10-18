## Use a simple multivariate normal for testing MCMC functionality

library(shinystan)
library(TMB)
compile('mvn.cpp')
dyn.load(dynlib('mvn'))

## devtools::install('C:/Users/Cole/shinystan')

d <- 2
covar <- diag(d) #matrix(c(1,-.95, -.95,1), nrow=2)
obj <- MakeADFun(data=list(d=d, covar=covar), parameters=list(X=rnorm(d)))

set.seed(4)

rm(theta.trajectory)
nuts <- run_mcmc(obj=obj, nsim=40, eps=.05, algorithm='NUTS', chains=1,
                covar=NULL, max_doubling=7)
xx <- (theta.trajectory)
plot(xx[,1], xx[,2], type='l')


sso <- with(nuts, as.shinystan(samples, burnin=warmup, max_treedepth=6,
                    sampler_params=sampler_params, algorithm='NUTS'))
launch_shinystan(sso)

## Run model in Stan
library(rstan)
data <- list(covar=covar, Npar=d, x=rep(0, len=d))
inits <- lapply(1:3, function(x) list(X=rnorm(d)))
fit <- stan(file='mvn.stan', data=data, iter=1000, init=inits, chains=3,
            control=list(metric='unit_e', max_treedepth=6))
sso.stan <- as.shinystan(fit)
launch_shinystan(sso.stan)


