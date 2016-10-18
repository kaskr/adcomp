## Use a simple multivariate normal for testing MCMC functionality

library(shinystan)
library(TMB)
compile('mvn.cpp')
dyn.load(dynlib('mvn'))

## devtools::install('C:/Users/Cole/shinystan')

d <- 2
chains <- 1
nsim <- 10
eps <- .6
td <- 1
covar <- diag(d) #matrix(c(1,-.95, -.95,1), nrow=2)
obj <- MakeADFun(data=list(d=d, covar=covar), parameters=list(X=rnorm(d)))

nuts <- run_mcmc(obj=obj, nsim=nsim, eps=eps, algorithm='NUTS', chains=chains,
                covar=NULL, max_doubling=td)
sso <- with(nuts, as.shinystan(samples, burnin=warmup, max_treedepth=td,
             sampler_params=sampler_params, algorithm='NUTS'))
launch_shinystan(sso)

## Run model in Stan
library(rstan)
data <- list(covar=covar, Npar=d, x=rep(0, len=d))
inits <- lapply(1:chains, function(x) list(X=rnorm(d)))
temp <- stan(file='mvn.stan', data=data, iter=5, init=inits, chains=1)
fit <- stan(fit=temp, data=data, iter=nsim, init=inits, chains=chains,
            control=list(metric='unit_e', stepsize=eps, adapt_engaged=FALSE,
                         max_treedepth=td))
sso.stan <- as.shinystan(fit)
launch_shinystan(sso.stan)


table(nuts$sampler_params[[1]][,4])

## Test single trajectories
set.seed(33)
theta <- theta0 <- rnorm(2)
r <- r0 <- rnorm(length(theta),0,1)
j <- 8
v <- 1
u <- Inf
rm(theta.trajectory)
info <- as.environment(list(n.calls=0))
u <- .sample.u(theta=theta, r=r, fn=fn)
xx <- .buildtree(theta=theta, r=r, u=u, v=v, j=j-1, eps=.3,
                 theta0=theta0, r0=r0, fn=fn, gr=gr, info=info)
info$n.calls
nrow(theta.trajectory)
plot(theta.trajectory[,1], theta.trajectory[,2], type='b', xlim=c(-3,3), ylim=c(-3,3))
