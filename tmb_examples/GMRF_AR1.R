#/*##########################################################################*/
#' GMRF-AR1 Spacetime Example
#' ==========================================================================
#'
#' by Mark R Payne
#' DTU-Aqua, Charlottenlund, Denmark
#' mpa@aqua.dtu.dk
#'
#' Thu Aug  7 08:30:07 2014
#'
#' Creates and solves a simple space-time example using TMB. The spatial
#' dimension is modelled as a GMRF, while the temporal is modelled as an
#' AR1 process
#' 
#  This work is subject to a Creative Commons "Attribution" "ShareALike" License.
#  You are largely free to do what you like with it, so long as you "attribute" 
#  me for my contribution. See the fine print at the end for exact details.
#
#  To do:
#
#  Notes:
#/*##########################################################################*/

# ========================================================================
# Initialise system
# ========================================================================
cat(sprintf("\n%s\n","GMRF-AR1 Spacetime Example"))
cat(sprintf("Analysis performed %s\n\n",date()))

#Configure markdown style, do house cleaning
rm(list = ls(all.names=TRUE));  graphics.off();
start.time <- proc.time()[3]; options(stringsAsFactors=FALSE)

#Helper functions, externals and libraries
log.msg <- function(fmt,...) {cat(sprintf(fmt,...));
                              flush.console();return(invisible(NULL))}
library(TMB)
library(MASS)
library(lattice)

# ========================================================================
#'## Settings
# ========================================================================
spat.grid.n <- 20    #Grid is n x n
temp.grid.n <- 12    #Temporal steps
obs.sd <- 0.1     #Observational noise
spat.cor.len <- 10   #Spatial decorrelation length in the GMRF
tsteps.to.drop <- c(6,8)  #Not included in the model fit
set.seed(123)

# ========================================================================
#'## Setup synthetic data
#'The simulated data is of the form z = A sin(t) + e
#'where A is a GMRF field of amplitudes, t is time and e is normal noise
# ========================================================================
#Create GMRF correlation matrix
grd <- expand.grid(x=seq(spat.grid.n),
                          y=seq(spat.grid.n))
grd.dis    <- as.matrix(dist(grd))
cor.mat    <- exp(-grd.dis/spat.cor.len)

#Create amplitude field
grd$amp.field <- mvrnorm(n=1,mu=rep(0,spat.grid.n^2),Sigma = cor.mat)

#Now add time dimension
dat <- expand.grid(space.idx=1:nrow(grd),t=1:temp.grid.n)
dat <- data.frame(grd[dat$space.idx,],dat)
dat$spacetime.idx <- 1:nrow(dat)
dat$truth <- dat$amp.field * sin((dat$t-1)/temp.grid.n*2*pi)

#Add observational noise
dat$z <- dat$truth + rnorm(nrow(dat),mean = 0,sd = obs.sd)

# ========================================================================
#'## Setup GMRF and fit
# ========================================================================
#Remove some of the time steps to demonstrate the ability of the
#method to handle gaps
dat.fit <- subset(dat,!(t %in% tsteps.to.drop))

#Setup GMRF
datin.l <- list(spacegrid=t(grd[,1:2]),
                obs=dat.fit$z,
                spacetime_idx=dat.fit$spacetime.idx)
par.init <- list(phi_trans=0,logdelta=0,logscale=0,logsigma=0,
                 eta=array(0,dim=c(spat.grid.n,spat.grid.n,temp.grid.n)))

#Compile
comp.status <- compile("GMRF_AR1.cpp")
if(comp.status!=0) stop("Compilation error.")
dyn.load(dynlib("GMRF_AR1"))

#Build objective function
objfn <- MakeADFun(data=datin.l,parameters=par.init,random="eta",DLL="GMRF_AR1")

#Fit
fit <- do.call(optim,objfn)

# ========================================================================
#'## Do plotting
# ========================================================================
#Extract values
res <- sdreport(objfn)
dat$fitted <- res$par.random
dat$err  <- dat$fitted-dat$truth

#Plot
zlim <- max(pretty(c(0,abs(dat$amp.field),abs(dat$fitted))))
zats <- seq(-zlim,zlim,by=0.1)
cols <- colorRampPalette(c("red","white","blue"))
sub.text <- paste("No data included in fit for time steps:",
                  paste(tsteps.to.drop,collapse=" "))

levelplot(truth ~ x + y | sprintf("Time step %02i",t),data=dat,
          at=zats,as.table=TRUE,
          col.regions=cols,
          main="Space-time truth")

plot(z~truth,dat,main="Truth + noise")

levelplot(fitted ~ x + y | sprintf("Time step %02i",t),data=dat,
          at=zats,as.table=TRUE,
          col.regions=cols,
          main="Space-time fit",sub=list(sub.text,cex=0.8))

levelplot(err ~ x + y | sprintf("Time step %02i",t),data=dat,
          at=zats,as.table=TRUE,
          col.regions=cols,          
          main="Fitted - truth",
          sub=list(sub.text,cex=0.8))

plot(fitted~truth,dat,main="Comparison - fitted vs truth")
abline(a=0,b=1,lwd=2,col="red")

#Is the precision of the missing time steps any different to the
#others?
boxplot(err ~ t, dat,xlab="Time step",ylab="Error",
        main="Comparison of errors across time steps",
        sub=sub.text,cex.sub=0.8)
abline(h=0,col="red")

# ========================================================================
# Complete
# ========================================================================
#Close files
if(grepl("pdf|png|wmf",names(dev.cur()))) {dmp <- dev.off()}
log.msg("\nAnalysis complete in %.1fs at %s.\n",proc.time()[3]-start.time,date())

#' -----------
#' <small>*This work by Mark R Payne is licensed under a  Creative Commons
#' Attribution-NonCommercial-ShareAlike 3.0 Unported License. 
#' For details, see http://creativecommons.org/licenses/by-nc-sa/3.0/deed.en_US
#' Basically, this means that you are free to "share" and "remix" for 
#' non-commerical purposes as you see fit, so long as you "attribute" me for my
#' contribution. Derivatives can be distributed under the same or 
#' similar license.*</small>
#'
#' <small>*This work comes with ABSOLUTELY NO WARRANTY or support.*</small>
#'
#' <small>*This work should also be considered as BEER-WARE. For details, see
#' http://en.wikipedia.org/wiki/Beerware*</small>
#' 
#' -----------
#
# Fin
