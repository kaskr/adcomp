
require(TMB)
compile('check_simulations.cpp')
dyn.load(dynlib('check_simulations'))
set.seed(1)

.results <- list()
checkSimulation <- function(obj) {
    g <- obj$gr(obj$par) ## Gradient at true parameter
    H <- obj$he(obj$par) ## Hessian at true parameter
    stopifnot(all( eigen(H)$value > 0 ))
    ## g ~ N(0, H)
    chi.sq <- g %*% solve(H) %*% t(g)
    p.value <- 1 - pchisq(chi.sq, df=length(obj$par))
    ans <- data.frame(distribution=obj$env$data$distr,
                      p.value=p.value, PASSED = p.value>.05 )
    .GlobalEnv$.results <- rbind(.GlobalEnv$.results, ans)
    print(ans)
    stopifnot(ans$PASSED)
    invisible(ans)
}

n <- 100000

######################################################################

data <- list(distr="norm", n=n)
parameters <- list(mu=0, sd=1)
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################

data <- list(distr="gamma", n=n)
parameters <- list(shape=3, scale=2)
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################

data <- list(distr="pois", n=n)
parameters <- list(lambda=5)
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################

data <- list(distr="compois", n=n)
parameters <- list(mode=3, nu=0.5)
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################

data <- list(distr="tweedie", n=n)
parameters <- list(mu=10, phi=10, p=1.5)
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################

data <- list(distr="compois2", n=n)
parameters <- list(mean=3, nu=0.1)
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################

data <- list(distr="nbinom", n=n)
parameters <- list(size=5, prob=.1)
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################

data <- list(distr="nbinom2", n=n)
parameters <- list(mu=5, var=10)
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################

data <- list(distr="exp", n=n)
parameters <- list(rate=1)
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################

data <- list(distr="beta", n=n)
parameters <- list(shape1=0.5, shape2=0.5)
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################

data <- list(distr="f", n=n)
parameters <- list(df1=5, df2=5)
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################

data <- list(distr="logis", n=n)
parameters <- list(location=0, scale=1)
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################

data <- list(distr="t", n=n)
parameters <- list(df=5)
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################

data <- list(distr="weibull", n=n)
parameters <- list(shape=1, scale=1)
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################
data <- list(distr="AR1", n=n)
parameters <- list(phi = .2)
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################

data <- list(distr="ARk", n=n)
parameters <- list(phi = c(.9, .05))
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################

data <- list(distr="MVNORM", n=1000)
parameters <- list(phi = c(.9))
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################
data <- list(distr="SEPARABLE", n=n)
parameters <- list(phi1 = .2, phi2 = c(.9, .05))
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################
data <- list(distr="GMRF", n=1000)
parameters <- list(delta = .1)
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################
data <- list(distr="SEPARABLE_NESTED", n=1000)
parameters <- list(phi1 = .2, phi2 = .9, delta = .1)
obj <- MakeADFun(data, parameters)
checkSimulation(obj)

######################################################################

cat("\n"); print(.results)
