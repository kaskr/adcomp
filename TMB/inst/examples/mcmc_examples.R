require(TMB)
## Draw samples with the RWM, HMC and NUTS algorithms and compare.

## Run the simple example, so that obj and opt are loaded into workspace
runExample("simple")

## Run two gradient based algorithms, using adative step size (eps) for
## each.
rwm <- mcmc(obj=obj, nsim=500, algorithm='RWM', params.init=opt$par,
            alpha=.08, diagnostic=TRUE)
hmc <- mcmc(obj=obj, nsim=500, algorithm='HMC', L=1, params.init=opt$par,
            diagnostic=TRUE)
nuts <- mcmc(obj=obj, nsim=500, algorithm='NUTS', params.init=opt$par,
             diagnostic=TRUE)

## Look at the adaptation of eps
with(hmc, plot(epsbar))
with(nuts, plot(epsbar))

## See how they compare via ACF
par(mfrow=c(3,4))
for(i in 1:4) acf(rwm$par[,i])
for(i in 1:4) acf(hmc$par[,i])
for(i in 1:4) acf(nuts$par[,i])

## also look at effective size per time
min(coda::effectiveSize(rwm$par))/rwm$time
min(coda::effectiveSize(hmc$par))/hmc$time
min(coda::effectiveSize(nuts$par))/nuts$time

rm(list=c('obj', 'opt', 'rwm', 'hmc', 'nuts'))
