## Use a simple multivariate normal for testing MCMC functionality

library(shinystan)
library(TMB)
compile('mvn.cpp')
dyn.load(dynlib('mvn'))

## devtools::install('C:/Users/Cole/shinystan')

covar <- matrix(c(1,-.9, -.9,1), nrow=2)
obj <- MakeADFun(data=list(d=2, covar=covar), parameters=list(X=c(1,-1)))

hmc <- run_mcmc(obj=obj, nsim=1000, algorithm='HMC', chains=2, L=5,
                covar=NULL)
str(hmc)
sso <- as.shinystan(hmc$samples, burnin=500,
                    sampler_params=hmc$sampler_params, algorithm='HMC')
launch_shinystan(sso)

nuts <- run_mcmc(obj=obj, nsim=1000, eps=.5, algorithm='NUTS', chains=2,
                covar=NULL)
str(nuts)
sso <- as.shinystan(nuts$samples, burnin=500, max_treedepth=8,
                    sampler_params=nuts$sampler_params, algorithm='NUTS')
launch_shinystan(sso)

