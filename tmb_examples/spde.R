# Leukemia example from Lindgren et al 2011, JRSS-B
# http://www.r-inla.org/examples/case-studies/lindgren-rue-and-lindstrom-rss-paper-2011
# Uses INLA to set up grid but use TMB to fit model

library(TMB)
compile("spde.cpp")
dyn.load(dynlib("spde"))

# Sets INLA mesh. Read R-INLA documentation about how to do this
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
data$spde <- (inla.spde2.matern(mesh, alpha=2)$param.inla)[c("M0","M1","M2")]	# Objects to be
n_s = nrow(data$spde$M0)

parameters <- list(beta=c(-5.0,0,0,0,0),log_tau=-2.0,log_kappa=2.5,log_alpha=-1,x=rep(0.0,n_s))

# Phase 1: Fit model non-spatial model first to get good starting values
obj <- MakeADFun(data,parameters,map=list(x=factor(rep(NA,n_s))),DLL="spde")
L=c(-7,-1,-1,-1,-1,-3.0,2.0,log(0.1))
U=c(-4,1,1,1,1,-1.0,3.0,log(10.0))
opt1 <- nlminb(obj$par,obj$fn,obj$gr,lower=L,upper=U)

# Phase 2: Include spatial part. Use starting values from phase 1
obj <- MakeADFun(data,parameters,random="x",DLL="spde")
opt <- nlminb(opt1$par,obj$fn,obj$gr,lower=L,upper=U)
