## From package 'geoR':
matern <- function (u, phi, kappa) 
{
    if (is.vector(u)) 
        names(u) <- NULL
    if (is.matrix(u)) 
        dimnames(u) <- list(NULL, NULL)
    uphi <- u/phi
    uphi <- ifelse(u > 0, (((2^(-(kappa - 1)))/ifelse(0, Inf, 
        gamma(kappa))) * (uphi^kappa) * besselK(x = uphi, nu = kappa)), 
        1)
    uphi[u > 600 * phi] <- 0
    return(uphi)
}

## Test:
n <- 100
set.seed(123)
x <- 10 * seq(0, 1, length=n)
D <- as.matrix(dist(x))
C <- matern(D, phi=1.3, kappa=2.7)
x <- t(chol(C)) %*% rnorm(n)
data <- list(x=x, D=D)
parameters <- list(phi=.5, kappa=.5)
map <- NULL

require(TMB)
compile('matern.cpp')
dyn.load(dynlib('matern'))

################################################################################

model <- MakeADFun(data, parameters, map=map)
system.time( fit <- nlminb(model$par, model$fn, model$gr) )
system.time( rep <- sdreport(model) )
print(rep)
