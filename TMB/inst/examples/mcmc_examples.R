\dontrun{
## Draw samples with the RWM, HMC and NUTS algorithms and compare.

## Run the simple example, so that obj and opt are loaded into workspace
runExample("simple")

## Run RWM and two gradient based algorithms, using adative step size (eps)
## for each. Start from the MLE.
rwm <- run_mcmc(obj=obj, nsim=500*8, algorithm='RWM', params.init=opt$par,
            alpha=.08, diagnostic=TRUE)
## Thin it to better approximate the gradient methods
rwm$par <- rwm$par[seq(1, nrow(rwm$par), by=8),]
hmc <- run_mcmc(obj=obj, nsim=500, algorithm='HMC', L=8, params.init=opt$par,
            diagnostic=TRUE, eps=0.1)
nuts <- run_mcmc(obj=obj, nsim=500, algorithm='NUTS', params.init=opt$par,
             diagnostic=TRUE, eps=0.1)

## See how they compare via ACF
par(mfrow=c(3,4))
for(i in 1:4) acf(rwm$par[,i])
for(i in 1:4) acf(hmc$par[,i])
for(i in 1:4) acf(nuts$par[,i])
}
