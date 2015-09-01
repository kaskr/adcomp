require(TMB)

## Draw samples with the HMC and NUTS algorithms and compare.

## Run the simple example, obj and opt are loaded into workspace
runExample("simple")
obj$env$beSilent()
## A helper function to get an approximate epsilon.
find.epsilon(theta=opt$par, fn=function(x) -obj$fn(x), gr=function(x) -obj$gr(x))
## Run two gradient based algorithms
system.time(simple.hmc <-
    run_mcmc(obj=obj, nsim=1000, algorithm='HMC', L=5, eps=.1, params.init=opt$par))
system.time(simple.nuts <-
    run_mcmc(obj=obj, nsim=1000, algorithm='NUTS', eps=.1, params.init=opt$par))
## See how they compare
par(mfrow=c(2,4))
for(i in 1:4) acf(simple.hmc[,i])
for(i in 1:4) acf(simple.nuts[,i])

rm(list=c('obj', 'opt', 'simple.hmc', 'simple.nuts'))
