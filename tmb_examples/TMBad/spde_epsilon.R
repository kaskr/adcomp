library(TMB)
compile("spde_epsilon.cpp", supernodal=FALSE, framework="TMBad") ## IMPORTANT: use 'supernodal=TRUE' if possible!
dyn.load(dynlib("spde_epsilon"))

## get cached objects - See 'spde_mesh.R'
##  'inla_mesh'
##  'inla_spde'
##  'Leuk'
load("spde_mesh.RData")

## Fixed effects part of model
X <- model.matrix( ~ 1 + sex + age + wbc + tpi, data = Leuk)

data <- list(time       = Leuk$time,
             notcens    = Leuk$cens,
             meshidxloc = inla_mesh$idx$loc - 1,
             X          = as.matrix(X))

## SPDE part: builds 3 components of Q (precision matrix)
data$spde <- inla_spde$param.inla[c("M0","M1","M2")]  # Encapsulation of 3 matrices
n_s <- nrow(data$spde$M0)                             # Number of points in mesh (including supporting points)

## Fit model
parameters <- list(beta      = c(-5.0,0,0,0,0),
                   log_tau   = -2.0,
                   log_kappa = 2.5,
                   log_omega = -1,
                   x         = rep(0.0, n_s),
                   epsilon   = numeric(0) )
obj <- MakeADFun(data, parameters, random="x", DLL="spde_epsilon", intern=!TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)

## ========= Epsilon method the normal way
## Add epsilon parameter
parameters$epsilon <- 0
p <- c(opt$par, epsilon=0)
system.time(obj <- MakeADFun(data, parameters, random="x", DLL="spde_epsilon"))
system.time(g <- obj$gr(p)) ## SLOW!

## ========= Epsilon method the new way
## detect sparse + lowrank hessian
control <- list(sparse=TRUE, lowrank=TRUE, trace=TRUE)
system.time(obj <- MakeADFun(data, parameters, random="x", DLL="spde_epsilon", intern=TRUE, inner.control=control))
system.time(g2 <- obj$gr(p)) ## FAST!

## Check error:
g[length(g)] - g2[length(g2)]
