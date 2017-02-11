#/*##########################################################################*/
#' GMRF with TMB Example
#' ==========================================================================
#'
#' by Mark R Payne
#' DTU-Aqua, Charlottenlund, Denmark
#' mpa@aqua.dtu.dk
#'
#' Tue Aug  5 15:38:45 2014
#'
#' Solves the classic FELSPINE problem that is used as an example in the mgcv
#' package (particularly in relation to soap bubble smoothing) using a GMRF. 
#' The problem is solved here using TMB. This example is a modified version of the
#' example from smooth.construct.so.smooth.spec() in the mgcv package.
#
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
cat(sprintf("\n%s\n","GMRF with TMB"))
cat(sprintf("Analysis performed %s\n\n",date()))

#Configure markdown style, do house cleaning
rm(list = ls(all.names=TRUE));  graphics.off();
start.time <- proc.time()[3]; options(stringsAsFactors=FALSE)

#Helper functions, externals and libraries
log.msg <- function(fmt,...) {cat(sprintf(fmt,...));
                              flush.console();return(invisible(NULL))}
library(mgcv)
library(TMB)

# ========================================================================
#'## Settings
# ========================================================================
#Observational noise
obs.sd <- 1

# ========================================================================
#'## Setup synthetic data
# ========================================================================
#Extract boundary
fsb <- fs.boundary()

#Create an underlying grid 
#Based on mgcv::fs.boundary() example
dx<-0.2;dy<-0.2    #Grid steps
id.fmt <- "%i/%i"
x.vals <- seq(-1,4,by=dx)
y.vals<-seq(-1,1,by=dy)
dat <- expand.grid(x=x.vals,y=y.vals)
tru.mat <- matrix(fs.test(dat$x,dat$y),length(x.vals),length(y.vals))

#Evaluate the function on it
dat$truth <- as.vector(tru.mat)
dat$x.idx <- as.numeric(factor(dat$x,x.vals))
dat$y.idx <- as.numeric(factor(dat$y,y.vals))
dat <- subset(dat,!is.na(truth))

## plot boundary with truth and data locations
par(mfrow=c(2,2))
image(x.vals,y.vals,tru.mat,col=heat.colors(100),xlab="x",ylab="y",
      main="Truth")
contour(x.vals,y.vals,tru.mat,levels=seq(-5,5,by=.25),add=TRUE)
lines(fsb$x,fsb$y);

#Add noise to the truth
dat$z <- dat$truth + rnorm(nrow(dat),sd=obs.sd) ## add noise
plot(z~truth,dat,main="Truth + noisy observations")
abline(a=0,b=1,lwd=2)

# ========================================================================
#'## Setup GMRF and fit
# ========================================================================
#Get grid details
grd.cells <- t(which(!is.na(tru.mat),arr.ind = TRUE))

#Associate observation to a grid cell
grd.ids <- sprintf("%i/%i",grd.cells[1,],grd.cells[2,])
dat$obs.id <- sprintf("%i/%i",dat$x.idx,dat$y.idx)
dat$obs.idx <- as.numeric(factor(dat$obs.id,levels=grd.ids))
  
#Setup GMRF
datin.l <- list(grid=grd.cells,
                obs=dat$z,n_obs=nrow(dat),
                obs_idx=dat$obs.idx)
par.init <- list(logdelta=0,logscale=0,logSigma=0,
                 eta=rep(0,ncol(grd.cells)))

#Compile
comp.status <- compile("GMRF.cpp")
if(comp.status!=0) stop("Compilation error.")
dyn.load(dynlib("GMRF"))

#Build objective function
objfn <- MakeADFun(data=datin.l,parameters=par.init,random="eta",DLL="GMRF")

#Fit
fit <- do.call(optim,objfn)

# ========================================================================
#'## Do plotting
# ========================================================================
#And plot
res <- sdreport(objfn)

dat$fitted <- res$par.random
plt.mat <- tru.mat
plt.mat[] <- NA
plt.mat[cbind(dat$x.idx,dat$y.idx)] <- dat$fitted
image(x.vals,y.vals,plt.mat,
      col=heat.colors(100),xlab="x",ylab="y",main="Fitted values")
points(y~x,dat,pch=3,cex=0.25)

plot(fitted~truth,dat,main="Comparison - fitted vs truth")
abline(a=0,b=1,lwd=2)

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
