# Anisotropic version of leukemia example from Lindgren et al 2011, JRSS-B
# http://www.r-inla.org/examples/case-studies/lindgren-rue-and-lindstrom-rss-paper-2011
# Uses INLA to set up grid but use TMB to fit model

library(TMB)
compile("spde_aniso.cpp")
dyn.load(dynlib("spde_aniso"))

# Sets up model and mesh in INLA
require(splancs)
require(rgl)
require(INLA)
require(lattice)
data(Leuk)
Leuk$id = 1:dim(Leuk)[1]
loc = cbind(Leuk$xcoord,Leuk$ycoord)
loc = cbind(Leuk$xcoord, Leuk$ycoord)
bnd1 = inla.nonconvex.hull(loc, convex=0.05)
bnd2 = inla.nonconvex.hull(loc, convex=0.25)
mesh = inla.mesh.2d(
                  ## Data locations to use as location seeds:
                  loc=cbind(Leuk$xcoord, Leuk$ycoord),
                  ## Encapsulate data region:
                  boundary=list(bnd1, bnd2),
                  ## Refined triangulation,
                  ## minimal angles >=26 degrees,
                  ## interior maximal edge lengths 0.05,
                  ## exterior maximal edge lengths 0.2,
                  ## don't add input points closer than 0.05:
                  min.angle=24,
                  max.edge=c(0.05, 0.2),
                  cutoff=0.005,
                  ## Set to >=0 for visual (no effect Windows):
                  plot.delay=0.5
                  )
     
# Fixed effects part of model
X = model.matrix(~1+sex+age+wbc+tpi,data=Leuk)
data <- list(time=Leuk$time,notcens=Leuk$cens,meshidxloc=mesh$idx$loc-1,X=as.matrix(X))

# SPDE part
spde = inla.spde2.matern(mesh, alpha=2)

# ---------- Start code that prepare objects for anisotropy. Should not be modified when you make your own example
Dset = 1:2
# Triangle info
TV = mesh$graph$tv       	 # Triangle to vertex indexing
V0 = mesh$loc[TV[,1],Dset]   # V = vertices for each triangle
V1 = mesh$loc[TV[,2],Dset]
V2 = mesh$loc[TV[,3],Dset]
E0 = V2 - V1                      # E = edge for each triangle
E1 = V0 - V2
E2 = V1 - V0  
# Calculate Areas
TmpFn = function(Vec1,Vec2) abs(det( rbind(Vec1,Vec2) ))
Tri_Area = rep(NA, nrow(E0))
for(i in 1:length(Tri_Area)) Tri_Area[i] = TmpFn( E0[i,],E1[i,] )/2   # T = area of each triangle
# ---------- End code that prepare objects for anisotropy. 


data$spde <- list(
      "n_s"=spde$n.spde,
      "n_tri"=nrow(TV),
      "Tri_Area"=Tri_Area,
      "E0"=E0,
      "E1"=E1,
      "E2"=E2,
      "TV"=TV-1,
      "G0"=spde$param.inla$M0,
      "G0_inv"=as(diag(1/diag(spde$param.inla$M0)),"dgTMatrix"))

parameters <- list(beta=c(-5.0,0,0,0,0),log_tau=-2.0,log_kappa=2.5,ln_H_input=c(0,0),log_alpha=-1,x=rep(0.0,data$spde$n_s))

# Fit model non-spatial model first to get good starting values
obj <- MakeADFun(data,parameters,random="x",DLL="spde_aniso")
L=c(-7,-1,-1,-1,-1,-3.0,2.0,c(-10,-10),log(0.1))
U=c(-4,1,1,1,1,-1.0,3.0,c(10,10),log(10.0))
opt <- nlminb(obj$par,obj$fn,obj$gr,lower=L,upper=U)

