

# load libraries
library(INLA)
library(TMB)

# Settings
Species = 'Simulated_counts'
  MeshType = c("Samples", "Refined")[1] # Samples: faster; Refined: slower and more accurate approximation to GRMF

# Loop across specise
  # Read data
  source("spatial_data.R")
  Data = read.csv( file="spatial_gompertz_simulated_data.csv", header=TRUE)

  # Length of data
  n_years = length(unique(Data$Year))
  n_stations = length(unique(Data$Site))
  n_data = n_stations*n_years
  x_stations = Data[match(unique(Data$Site),Data$Site),'Lon..DDD.DDDDD.']
  y_stations = Data[match(unique(Data$Site),Data$Site),'Lat..DD.DDDDD.']

  # reformat data
  Ymat = matrix(NA, nrow=n_stations, ncol=n_years)
  for(YearI in 1:n_years){
  for(SiteI in 1:n_stations){
    Which = which(Data$Year==unique(Data$Year)[YearI] & Data$Site==unique(Data$Site)[SiteI])
    if(length(Which)>=1) Ymat[SiteI,YearI] = Data[Which[1],Species]
  }}
  # vectorize data
  Y = as.vector(Ymat)		  
  X = cbind( rep(1,n_data) )
  Site = as.vector(row(Ymat))
  Year = as.vector(col(Ymat))
  NAind = as.integer(ifelse(is.na(Y),1,0))

  # Build SPDE object using INLA
  if(MeshType=="Samples") mesh = inla.mesh.create( cbind(x_stations, y_stations), plot.delay=NULL, extend=list(n=8,offset=-0.15), refine=F )  # loc_samp
  if(MeshType=="Refined") mesh = inla.mesh.create( cbind(x_stations, y_stations), plot.delay=NULL, extend=list(n=8,offset=-0.15), refine=list(min.angle=26) )  # loc_samp  ;  ,max.edge.data=0.08,max.edge.extra=0.2
  spde = inla.spde2.matern(mesh,alpha=2)

  # Settings
  newtonOption(smartsearch=TRUE)
  setwd( TmbFile )

  # Run spatial model
  dyn.load( dynlib("spatial_gompertz") )
    Data = list(n_data=n_stations*n_years, Y=Y, NAind=NAind, n_stations=n_stations, meshidxloc=mesh$idx$loc-1, n_years=n_years, n_p=ncol(X), X=X, G0=spde$param.inla$M0, G1=spde$param.inla$M1, G2=spde$param.inla$M2)
    Parameters = list(alpha=c(0.0), phi=0.0, log_tau_E=0.0, log_tau_O=0.0, log_kappa=0.0,	rho=0.5, Epsilon_input=matrix(rnorm(spde$n.spde*n_years),nrow=spde$n.spde,ncol=n_years), Omega_input=rnorm(spde$n.spde))
    Random = c("Epsilon_input","Omega_input")
    obj <- MakeADFun(data=Data, parameters=Parameters, random=Random, hessian=FALSE)
    obj$fn(obj$par)

    # Run optimizer
    obj$control <- list(trace=1,parscale=rep(1,13),REPORT=1,reltol=1e-12,maxit=100)
    opt = nlminb(obj$par, obj$fn, obj$gr, lower=c(rep(-20,2),rep(-10,3),-0.999), upper=c(rep(20,2),rep(10,3),0.999), control=list(eval.max=1e4, iter.max=1e4))
    SD = try( sdreport(obj) )
    Report = obj$report()

  # Range of correlation (Lindgren and Rue 2013, immediately before Eq. 4)
    Nu = 1
    sqrt(8*Nu)/exp(opt$par['log_kappa'])
  # Marginal Standard Deviation  (Lindgren and Rue 2013 Eq. 2)
    1 / sqrt(4*pi*exp(2*opt$par['log_tau_E'])*exp(2*opt$par['log_kappa']))

