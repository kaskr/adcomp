## randomwalkvalidation
##
## Estimate and validate a random walk model with and without drift
##
## Compare Thygesen et al (submitted, 2016): Validation of state space models
## fitted as mixed effects models
##
## Uffe Høgsbro Thygesen and Kasper Kristensen, 2016

library(TMB)
compile("randomwalkvalidation.cpp")
dyn.load(dynlib("randomwalkvalidation"))

PlotToFile <- FALSE

## For paper versions of graphs
figwidth <- 5
figheight <- 5

## For reproducible results
set.seed(123)

## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

## Matlab-style stem plot
stem <- function(x,y,pch=16,linecol=1,clinecol=1,...){
if (missing(y)){
    y = x
    x = 1:length(x) }
    plot(x,y,pch=pch,...)
    for (i in 1:length(x)){
       lines(c(x[i],x[i]), c(0,y[i]),col=linecol)
    }
    lines(c(x[1]-2,x[length(x)]+2), c(0,0),col=clinecol)
}

## This function plots a time series and to the side (left or right) a
## histogram
timeseries.w.histogram <- function(x,breaks=-7:7,histleft=TRUE,col="black",...)
{
        if(max(abs(x),na.rm=TRUE)>max(abs(breaks)))
        {
            breaks <- NULL
        }
        
        if(is.null(breaks))
            {
                maxAbsX <- ceiling(max(abs(x),na.rm=TRUE)) 
                breaks <- seq(-maxAbsX,maxAbsX)
           }

        if(length(breaks)==1)
            {
                breaks <- seq(-breaks,breaks)                
            }

        plot(x,ylim=range(breaks),col=col,...)
        h <- hist(x[is.finite(x)],plot=FALSE,breaks=breaks)
        xscale <- 0.25*length(x)/max(h$density)

        xcoor <- c(0,h$density*xscale)
        if(!histleft)
            {
                xcoor <- length(x)-xcoor
            }
        
        y <- h$breaks
        ny <- length(y)
        y <- c(y[ny],rep(y[-ny],rep(2,ny-1)),y[ny])

        col = add.alpha(col,0.2)
        
        polygon(rep(xcoor,rep(2,length(xcoor))),
                y,border=NA, col=col,... )
    }

## Simulate data with these parameters
mu <- 0.75
sigma <- 1
s <- 1
huge <- 1e3

## Simulate random track
nT <- 100
X <- c(0,cumsum(rnorm(nT-1,mean=mu,sd=sigma)))

## Simulate random measurements
Y <- X + rnorm(nT,sd=s)

data <- list(y=Y,huge=huge)
parameters <- list(x=X,mu=0,logsigma=log(sigma),logs=log(s))

## Estimate states and parameters under H0: mu=0
obj0 <- MakeADFun(data,parameters,random=c("x"),DLL="randomwalkvalidation",map=list(mu=factor(NA)))
opt0 <- do.call("optim",obj0)
sdr0 <- sdreport(obj0)
estX0 <- summary(sdr0,"random")

## Estimate states and parameters under H1: mu != 0
obj1 <- MakeADFun(data,parameters,random=c("x"),DLL="randomwalkvalidation")
opt1 <- do.call("optim",obj1)
sdr1 <- sdreport(obj1)
estX1 <- summary(sdr1,"random")

## Plot true states, estimated states, measured states
if(PlotToFile) pdf(file="X.pdf",width=figwidth,height=figheight) else dev.new()

plot(X,type="l",xlab="Time",ylab="X",lwd=3)
points(Y)
lines(estX0[,1],lwd=1,col="red")
lines(estX1[,1],lwd=1,col="green")
legend("topleft",
       legend=c("True state X","Estimated state X (H0)",
           "Estimated state X (H1)","Measured state Y"),
       lty=c("solid","solid","solid",NA),
       lwd=c(3,1,1,NA),
       pch=c(NA,NA,NA,"o"),
       col=c("black","red","green","black")
       )
par(fig=c(0.5,0.98,0.08,0.7),new=TRUE)
plot(X,type="l",xlab=NA,ylab=NA,xlim=c(30,34),ylim=range(X[30:35]),lwd=3)
points(Y)
lines(estX0[,1],lwd=1,col="red")
lines(estX1[,1],lwd=1,col="green")

if(PlotToFile) dev.off()

### Compute normalized smoothing residuals
### Note: the point is that these are unsuitable for validation

## Plot a.c.f. of normalized smoothing residuals
resid0 <- Y-estX0[,1]
resid1 <- Y-estX1[,1]

Norm.resid0 <- resid0 / estX0[,2]
Norm.resid1 <- resid1 / estX1[,2]

acf0 <- acf(Norm.resid0,plot=FALSE)
acf1 <- acf(Norm.resid1,plot=FALSE)

if(PlotToFile) pdf(file="Acf.pdf",width=figwidth,height=figheight) else dev.new()

xlim <- range(c(acf0$lag,acf1$lag))
ylim <- range(c(acf0$acf,acf1$acf))
dx <- 0.1
stem(acf0$lag-dx,acf0$acf,col="red",linecol="red",
     ylim=ylim,xlim=xlim,xlab="Lag",ylab="A.c.f.")
par(new=TRUE)
stem(acf1$lag+dx,acf1$acf,col="green",linecol="green",
     ylim=ylim,xlim=xlim,xlab=NA,ylab=NA)
ci <- 0.95
clim <- qnorm((1 + ci)/2)/sqrt(acf0$n.used)
lines(xlim,rep(clim,2),lty="dashed",col="grey")
lines(xlim,-rep(clim,2),lty="dashed",col="grey")

if(PlotToFile) dev.off()

## Plot smoothing residuals 
if(PlotToFile) pdf(file="Resid.pdf",width=figwidth,height=figheight) else dev.new()
ylim <- max(abs(c(resid0,resid1)))
breaks <- pretty(c(-1,1)*ylim,n=20)
                 
timeseries.w.histogram(resid0,xlab="Time",ylab="Naive residuals",
                       histleft=TRUE,breaks=breaks,col=rgb(1,0,0))
par(new=TRUE)
timeseries.w.histogram(resid1,xlab=NA,ylab=NA,
                       histleft=FALSE,breaks=breaks,col=rgb(0,1,0))
abline(h=0)
if(PlotToFile) dev.off()

print(paste("Mean resid0 = ",mean(resid0)))
print(paste("Sdev resid0 = ",sqrt(var(resid0))))

print(paste("Mean resid0 = ",mean(resid1)))
print(paste("Sdev resid0 = ",sqrt(var(resid1))))

### Validation using one sample from the posterior
require(MASS)

C0 <- solve(obj0$env$spHess(random=TRUE))  ## Covariance matrix of random effects
Xr0 <- mvrnorm(1,estX0[,1],C0)             ## Generate one sample of random effects
W0 <- diff(Xr0)/exp(opt0$par["logsigma"])  ## Compute corresponding driving process noise

C1 <- solve(obj1$env$spHess(random=TRUE))  ## .. repeat for H1
Xr1 <- mvrnorm(1,estX1[,1],C1)
W1 <- (diff(Xr1)-opt1$par["mu"])/exp(opt1$par["logsigma"])

## Plot process noise as time series
if(PlotToFile) pdf(file="Sampled-increments.pdf",width=figwidth,height=figheight) else dev.new()

breaks <- pretty(range(c(W0,W1)),15)

timeseries.w.histogram(W0,xlab="Time",ylab="Normalized process noise",
                       histleft=TRUE,breaks=breaks,col=rgb(1,0,0))
par(new=TRUE)
timeseries.w.histogram(W1,xlab=NA,ylab=NA,
                       histleft=FALSE,breaks=breaks,col=rgb(0,1,0))
abline(h=0)
if(PlotToFile) dev.off()

### One-step predictions

## Generate one step predictions with the models fitted under H0 and H1
predict0  <- oneStepPredict(obj0,observation.name="y",method="fullGaussian")
predict1  <- oneStepPredict(obj1,observation.name="y",method="fullGaussian")

## Plot both prediction residuals in one plote
if(PlotToFile)  pdf(file="Pred-resid.pdf",width=figwidth,height=figheight) else dev.new()

breaks <- pretty(range(predict0$residual,predict1$residual,na.rm=TRUE),15)

timeseries.w.histogram(predict0$residual,xlab="Time",
                       ylab="Normalized prediction errors",
                       histleft=TRUE,breaks=breaks,col=rgb(1,0,0))

par(new=TRUE)
timeseries.w.histogram(predict1$residual,xlab=NA,ylab=NA,
                       histleft=FALSE,breaks=breaks,col=rgb(0,1,0))
abline(h=0)
if(PlotToFile) dev.off()

## Test if the prediction errors are unbiased
print(TestPred1 <- anova(lm(predict1$residual ~ 0),lm(predict1$residual ~ 1)))
print(TestPred0 <- anova(lm(predict0$residual ~ 0),lm(predict0$residual ~ 1)))

