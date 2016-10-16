## Use a simple multivariate normal for testing MCMC functionality

library(shinystan)
library(TMB)
compile('mvn.cpp')
dyn.load(dynlib('mvn'))

## devtools::install('C:/Users/Cole/shinystan')

covar <- matrix(c(1,-.95, -.95,1), nrow=2)
obj <- MakeADFun(data=list(d=2, covar=covar), parameters=list(X=c(1,-1)))

hmc <- run_mcmc(obj=obj, nsim=1000, algorithm='HMC', chains=3, L=100,
                covar=NULL, eps=NULL)
sso <- as.shinystan(hmc$samples, burnin=500,
                    sampler_params=hmc$sampler_params, algorithm='HMC')
launch_shinystan(sso)

set.seed(4)
nuts <- run_mcmc(obj=obj, nsim=1000, eps=NULL, algorithm='NUTS', chains=2,
                covar=covar, max_doubling=8)
str(nuts)
sso <- as.shinystan(nuts$samples, burnin=500, max_treedepth=8,
                    sampler_params=nuts$sampler_params, algorithm='NUTS')
launch_shinystan(sso)

## Run model in Stan
library(rstan)
data <- list(covar=covar, Npar=2, x=rep(0, len=2))
inits <- list(list(X=c(0,0)),list(X=c(0,0)),list(X=c(0,0)))
fit <- stan(file='mvn.stan', data=data, iter=2000, init=inits, chains=3,
            control=list(metric='unit_e'))
sso.stan <- as.shinystan(fit)
launch_shinystan(sso)


