## Use a simple multivariate normal for testing MCMC functionality

library(shinystan)
library(TMB)
library(rstan)
compile('mvn.cpp')
dyn.load(dynlib('mvn'))

## devtools::install('C:/Users/Cole/shinystan')

d <- 10
chains <- 5
nsim <- 2000
eps <- 1
td <- 8
covar <- diag(d) #matrix(c(1,-.95, -.95,1), nrow=2)
#diag(covar) <- rexp(d, 2)
obj <- MakeADFun(data=list(d=d, covar=covar), parameters=list(X=rnorm(d)))

nuts <- run_mcmc(obj=obj, nsim=nsim, eps=NULL, algorithm='NUTS', chains=chains,
                covar=NULL, max_doubling=td)
sso <- with(nuts, as.shinystan(samples, burnin=warmup, max_treedepth=td,
             sampler_params=sampler_params, algorithm='NUTS'))
launch_shinystan(sso)
temp <- as.data.frame(sso@sampler_params[[1]])
all(temp$n_leapfrog__ <= 2^temp$treedepth__ -1)
all(temp$n_leapfrog__ >= 2^(temp$treedepth__-1) -1)
table(temp[,3])/nrow(temp)
plot(2^temp[,3]-1, temp[,4]); abline(a=0, b=1)
plot(nuts$sampler_params[[1]][,1])
## plot(nuts$sampler_params[[1]][,3])

## Run model in Stan
data <- list(covar=covar, Npar=d, x=rep(0, len=d))
inits <- lapply(1:chains, function(x) list(X=rnorm(d)))
temp <- stan(file='mvn.stan', data=data, iter=5, init=inits, chains=1)
fit <- stan(fit=temp, data=data, iter=nsim, init=inits, chains=chains,
            control=list(metric='unit_e', stepsize=eps, adapt_engaged=TRUE,
                         max_treedepth=td))
sso.stan <- as.shinystan(fit)
temp2 <- sso.stan@sampler_params[[1]]
table(temp2[,3])/nrow(temp2)
plot(2^temp2[,3]-1, temp2[,4]); abline(a=0, b=1)
launch_shinystan(sso.stan)
all(temp$n_leapfrog__ <= 2^temp$treedepth__ -1)
all(temp$n_leapfrog__ >= 2^(temp$treedepth__-1) -1)

lapply(sso@sampler_params, function(x) table(x[,3]))
lapply(sso.stan@sampler_params, function(x) table(x[,3]))

## Test single trajectories
set.seed(3)
theta <- theta0 <- rnorm(d)
r <- r0 <- rnorm(length(theta),0,1)
j <- 8
v <- 1
u <- Inf


rm(theta.trajectory)
u <- .sample.u(theta=theta, r=r, fn=fn)
info <- as.environment(list(n.calls=0))
xx <- .buildtree(theta=theta, r=r, u=u, v=v, j=12, eps=.1,
                 theta0=theta0, r0=r0, fn=fn, gr=gr, info=info)
info$n.calls
plot(theta.trajectory[,2], theta.trajectory[,3], type='b', xlim=c(-3,3), ylim=c(-3,3))
points(theta[1], theta[2], pch=16)
with(xx, points(theta.prime[1], theta.prime[2], pch=16, col='red', cex=.5))

with(xx, (theta.plus-theta.minus) %*% r.minus)
drop(with(xx, crossprod(theta.plus-theta.minus, r.minus)))
cbind(step=0, t(c(1,1)))

xx <- t(sapply(1:5000, function(i) .buildtree(theta=theta, r=r, u=u, v=v, j=j-1, eps=.3,
                 theta0=theta0, r0=r0, fn=fn, gr=gr, info=info)$theta.prime))
barplot(table(xx[,1]))
