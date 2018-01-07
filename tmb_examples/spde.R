# Illustration SPDE/INLA approach to spatial modelling via Matern correlation function
# Leukemia example from Lindgren et al 2011, JRSS-B
# http://www.r-inla.org/examples/case-studies/lindgren-rue-and-lindstrom-rss-paper-2011

library(TMB)
compile("spde.cpp")
dyn.load(dynlib("spde"))

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

parameters <- list(beta      = c(-5.0,0,0,0,0),
                   log_tau   = -2.0,
                   log_kappa = 2.5,
                   log_omega = -1,
                   x         = rep(0.0, n_s) )

## Phase 1: Fit non-spatial part first to get good starting values for fixed effects
not_phase1 <- list(log_tau   = as.factor(NA),
                   log_kappa = as.factor(NA),
                   x         = factor(rep(NA, n_s)) )
obj <- MakeADFun(data, parameters, map=not_phase1, DLL="spde")
opt1 <- nlminb(obj$par, obj$fn, obj$gr)

## Modify starting values after phase 1
parameters <- list(beta      = opt1$par[1:5],
                   log_tau   = -2.0,
                   log_kappa =  2.5,
                   log_omega = opt1$par["log_omega"],
                   x         = rep(0.0, n_s))

## Phase 2: Include spatial part. Use starting values from phase 1
obj <- MakeADFun(data, parameters, random="x", DLL="spde")
L   <- c(-7, -1, -1, -1, -1, -3.0, 2.0, log(0.1) )
U   <- c(-4,  1,  1,  1,  1, -1.0, 3.0, log(10.0))
opt <- nlminb(obj$par, obj$fn, obj$gr, lower=L, upper=U)

# Calculate standard deviations, and extract rho
Rep <- sdreport(obj)
rho_est <- summary(Rep,"report")
