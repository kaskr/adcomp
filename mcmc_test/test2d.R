require(TMB)

compile('test2d.cpp')
dyn.load(dynlib('test2d'))

## Image of density
plotSurface <- function(model){
    x <- seq(-5,5,length=51)
    e <- expand.grid(x=x,y=x)
    f <- apply(e,1,model$f)
    dim(f) <- rep(length(x),2)
    image(x,x,-exp(-f))
}

## Get independent sample
getOneSample <- function(method, plot = FALSE){
    sim <- switch(method,
                  "HMC"  = mcmc(model, 1000,  method, L=5, eps=.5),
                  "NUTS" = mcmc(model, 100,   method, eps=.5),
                  "RWM"  = mcmc(model, 10000, method) )
    ans <- tail(sim, 1)
    if(plot) points(ans)
    unlist(ans)
}

## Simulate independent draws from target:
set.seed(123)
nsim <- 100

## Tests all the mcmc algorithms
ans <- list()
for(example in c(0, 1)){
    parameters <- list(x=0, y=0)
    model <- MakeADFun(list(example=example), parameters)
    for(method in c("HMC", "NUTS", "RWM")){
        dev.new(); layout(t(1:2))
        plotSurface(model)
        system.time( sim <- replicate(nsim, getOneSample(method, TRUE)) )
        ## Transform the draws to iid uniform variates and test for
        ## correctness:
        uv <- apply(sim, 2, function(xy) model$report(xy)$uv)
        plot(seq(0, 1, length = 2*nsim), sort(uv)); abline(0,1)
        title(paste0(method,example))
        ans[[paste0(method,example)]] <- list(
            uv,
            cor(t(uv)),
            ks.test(uv, "punif")
            )
    }
}
