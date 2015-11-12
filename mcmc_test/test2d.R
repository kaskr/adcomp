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
                  "HMC"  = mcmc(model, 1000,  method, L=5, diagnostic=TRUE),
                  "NUTS" = mcmc(model, 100,   method, diagnostic=TRUE),
                  "RWM"  = mcmc(model, 10000, method, diagnostic=TRUE) )
    ans <- tail(sim$par,1)
    if(plot) points(ans)
    eps <- if(method=="RWM") 0 else tail(sim$epsbar,1)
    data.frame(ans, eps=eps)
}

## Simulate independent draws from target:
set.seed(123)
nsim <- 100

## Tests all the mcmc algorithms
eps <- ans <- list()
for(example in c(0, 1)){
    dev.new(); par(mfrow=c(3,2))
    parameters <- list(x=0, y=0)
    model <- MakeADFun(list(example=example), parameters)
    for(method in c("HMC", "NUTS", "RWM")){
        plotSurface(model)
        system.time( sim <- lapply(1:nsim, function(x) getOneSample(method, TRUE)) )
        ## Transform the draws to iid uniform variates and test for
        ## correctness:
        sim <- do.call(rbind, sim)
        uv <- apply(sim[,1:2], 1, function(xy) model$report(xy)$uv)
        plot(seq(0, 1, length = 2*nsim), sort(uv)); abline(0,1)
        title(paste0(method,example))
        ans[[paste0(method,example)]] <- list(
            uv,
            cor(t(uv)),
            ks.test(uv, "punif")
            )
        eps[[paste0(method,example)]] <- sim$eps
    }
}
dev.new()
## Make sure the adaptive eps algorithms are converging consistently
eps <- data.frame(do.call(cbind, eps))[,c(-3,-6)]
boxplot(eps)
apply(eps, 2, function(x) sum(!is.finite(x)))
