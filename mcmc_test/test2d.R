require(TMB)

compile('test2d.cpp')
dyn.load(dynlib('test2d'))

## Image of density
plotSurface <- function(model){
    x <- seq(-5,5,length=51)
    e <- expand.grid(x=x,y=x)
    f <- apply(e,1,model$f)
    dim(f) <- rep(length(x),2)
    image(x,x,-exp(-f), axes=FALSE, ann=FALSE); box()
}

## Get independent sample
getOneSample <- function(method, plot = FALSE, covar=NULL){
    sim <-
     switch(method,
       "HMC" =run_mcmc(model, 1000,  method, L=5, diagnostic=TRUE, covar=covar),
       "NUTS"=run_mcmc(model, 100,   method, diagnostic=TRUE, covar=covar),
       "RWM" =run_mcmc(model, 10000, method, diagnostic=TRUE, covar=covar))
    ans <- tail(sim$par,1)
    if(plot) points(ans)
    eps <- if(method=="RWM") 0 else tail(sim$epsbar,1)
    data.frame(ans, eps=eps)
}

## Simulate independent draws from target:
set.seed(123)
nsim <- 200

## xx <- run_mcmc(model, nsim=200, algorithm='HMC', L=5, params.init=c(5,5), covar=covar)


## Tests all the mcmc algorithms with a specified covariance
covar <- matrix(c(1, .5, .5, 1), nrow=2)
eps <- ans <- list()
for(example in c(0, 1)){
    dev.new(); par(mfrow=c(6,2), mar=.1*c(1,1,1,1))
    for(co in list(NULL, covar)){
        cov.label <- ifelse(is.null(co),"nocovar","covar")
        parameters <- list(x=0, y=0)
        model <- MakeADFun(list(example=example), parameters)
        for(method in c("HMC", "NUTS", "RWM")){
            scenario <- paste0(method,example,cov.label)
            plotSurface(model); title(scenario, line=-1)
            system.time( sim <- lapply(1:nsim, function(x)
                getOneSample(method, TRUE, covar=co)) )
            ## Transform the draws to iid uniform variates and test for
            ## correctness:
            sim <- do.call(rbind, sim)
            uv <- apply(sim[,1:2], 1, function(xy) model$report(xy)$uv)
            plot(seq(0, 1, length = 2*nsim), sort(uv), axes=FALSE,
                 ann=FALSE); box()
            abline(0,1); title(scenario, line=-1)
            ans[[scenario]] <- list(
                uv,
                cor(t(uv)),
                ks.test(uv, "punif")
                )
            eps[[scenario]] <- sim$eps
        }
    }
}

## Make sure the adaptive eps algorithms are converging consistently
dev.new()
boxplot(eps)
sapply(eps, function(x) sum(!is.finite(x)))
