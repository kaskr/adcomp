# Inference in a linear scalar stochastic differential equation
#
# dX = - lambda*X*dt + sigmaX*dB
#
# from discrete observations
#
# Y(i) = X(t(i)) + e(i)
#
# where e(i) is N(0,sigmaY^2)
#
# Latent variables are the states X
#
# We use Euler approximation to evalaute transition densities. The time mesh for this
# discretization is finer than the sample interval, i.e. some (many) states are unobserved.
#
# The code first simulates a data set, then reestimates parameters and states

set.seed(1);  # For reproducible results
library(TMB)

lambda <- -1   # Rate parameter in the SDE				     
sigmaX <- 0.1  # log(sigmaX) where sigmaX is noise intensity in the SDE   
sigmaY <- 1e-2 # log(sigmaY) where sigmaY is std.dev. on measurement error

x0 <- 0.5       # Initial state

f <- function(x) lambda*x
g <- function(x) sigmaX

Tsim <- 0.1     # Time step for Euler scheme
T <- 50         # Duration of simulation
Tobs <- 1       # Sample time interval. Should be divisible with Tsim

# "Generic" function to sample a sample path of an SDE using Euler
euler <- function(x,f,g,tvec,dB=NULL){
  X <- numeric(length(tvec))
  X[1] <- x0
  dt <- diff(tvec)
  if(is.null(dB)) dB <- rnorm(length(dt),sd=sqrt(dt))
  for(i in 1:(length(tvec)-1))
    X[i+1] <- X[i] + f(X[i])*dt[i] + g(X[i])*dB[i]
  return(X)
}

# Simulate trajectory
tsim <- seq(0,T,Tsim)
Xsim <- euler(x0,f,g,tsim)

# Measurements with sample interval Tobs
iobs <- seq(1,length(tsim),round(Tobs/Tsim))

# Generate random measurements
Y <- rnorm(length(iobs),mean=(Xsim[iobs]),sd = sigmaY)

dyn.load(dynlib("sde_linear"))

# Data for TMB
data <- list(tsim=tsim,iobs=iobs-1,Y=Y)

# Initial geuess on latent variables
X <- euler(x0,f,g,tsim,dB=numeric(length(tsim)-1))
dB <- numeric(length(tsim)-1)

parameters <- list(
                   X=X,
                   lambda=lambda,
                   logsX=log(sigmaX),
                   logsY=log(sigmaY)
                   )

# The inner problem is quadratic
newtonOption(smartsearch=FALSE)

obj <- MakeADFun(data,parameters,random=c("X"),DLL="sde_linear")

# Estimate parameters and latent variables
system.time(opt <- nlminb(obj$par,obj$fn,obj$gr))
sdr <- sdreport(obj)

# Get predictions of states with std.dev.
pred <- summary(sdr,"random")
Xpred <- pred[,1]
Xsd <- pred[,2]

# Setup plot
plot(tsim,Xsim,type="n",xlab="Time t",ylab="State x")

# Confidence region for states
polygon(c(tsim,rev(tsim)),c(Xpred+1.96*Xsd,rev(Xpred-1.96*Xsd)),col="grey",border=NA)

# Plot true path
lines(tsim,Xsim,type="l")

# Add predictions
lines(tsim,Xpred,pch=16,col="red")

# Add measurements
points(tsim[iobs],Y)

legend("topright",legend=c("True","Measured","Smoothed"),lty=c("solid",NA,"solid"),pch=c(NA,1,NA),col=c("black","black","red"))
