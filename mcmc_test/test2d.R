require(TMB)
parameters <- list(x=0, y=0)
compile('test2d.cpp')
dyn.load(dynlib('test2d'))
model <- MakeADFun(list(example=0), parameters)

## Image of density
x <- seq(-5,5,length=51)
e <- expand.grid(x=x,y=x)
f <- apply(e,1,model$f)
dim(f) <- rep(length(x),2)
image(x,x,f)

## FIXME: Also test the other mcmc algorithms
getOneSample <- function(plot = FALSE){
    sim <- run_mcmc(model, 10000, "HMC", L=5, eps=.1)
    ans <- tail(sim,1)
    if(plot) points(ans)
    unlist(ans)
}

## Simulate independent draws from target:
set.seed(123)
nsim <- 100
system.time( sim <- replicate(nsim, getOneSample(TRUE)) )

## Transform the draws to iid uniform variates and test for
## correctness:
uv <- apply(sim, 2, function(xy) model$report(xy)$uv)
plot(seq(0, 1, length = 2*nsim), sort(uv)); abline(0,1)
cor(t(uv))
ks.test(uv, "punif")
