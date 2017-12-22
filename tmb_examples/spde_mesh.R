## Sets INLA mesh. Please refer to R-INLA documentation
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
inla_mesh = inla.mesh.2d(
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

# INLA SPDE object
inla_spde = inla.spde2.matern(inla_mesh, alpha=2)

## save(inla_mesh, inla_spde, Leuk, file="spde_mesh.RData")
