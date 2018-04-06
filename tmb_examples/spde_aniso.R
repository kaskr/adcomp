# Anisotropic version of leukemia example from Lindgren et al 2011, JRSS-B
# http://www.r-inla.org/examples/case-studies/lindgren-rue-and-lindstrom-rss-paper-2011
# Uses INLA to set up grid but use TMB to fit model

library(TMB)
library(Matrix)
compile("spde_aniso.cpp")
dyn.load(dynlib("spde_aniso"))

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

# ---------- Start code that prepare objects for anisotropy. Should not be modified when you make your own example
Dset = 1:2
# Triangle info
TV = inla_mesh$graph$tv           # Triangle to vertex indexing
V0 = inla_mesh$loc[TV[,1],Dset]   # V = vertices for each triangle
V1 = inla_mesh$loc[TV[,2],Dset]
V2 = inla_mesh$loc[TV[,3],Dset]
E0 = V2 - V1                      # E = edge for each triangle
E1 = V0 - V2
E2 = V1 - V0  
# Calculate Areas
TmpFn = function(Vec1, Vec2) abs(det( rbind(Vec1, Vec2) ))
Tri_Area = rep(NA, nrow(E0))
for(i in 1:length(Tri_Area)) Tri_Area[i] = TmpFn( E0[i,],E1[i,] )/2   # T = area of each triangle
# ---------- End code that prepare objects for anisotropy. 

data$spde <- list(
    "n_s"      = inla_spde$n.spde,
    "n_tri"    = nrow(TV),
    "Tri_Area" = Tri_Area,
    "E0"       = E0,
    "E1"       = E1,
    "E2"       = E2,
    "TV"       = TV - 1,
    "G0"       = inla_spde$param.inla$M0,
    "G0_inv"   = as(diag(1/diag(inla_spde$param.inla$M0)), "dgTMatrix"))

parameters <- list(beta       = c(-5.0,0,0,0,0),
                   log_tau    = -2.0,
                   log_kappa  =  2.5,
                   ln_H_input = c(0,0),
                   log_alpha  = -1,
                   x          = rep(0.0, data$spde$n_s))

## Fit model
obj <- MakeADFun(data, parameters, random="x", DLL="spde_aniso")
L <- c(-7, -1, -1, -1, -1, -3.0, 2.0, c(-10,-10), log(0.1) )
U <- c(-4,  1,  1,  1,  1, -1.0, 3.0, c( 10, 10), log(10.0))
opt <- nlminb(obj$par, obj$fn, obj$gr, lower=L, upper=U)

## Report
sdr <- sdreport(obj)
