library(TMB)

## compile("hmm.cpp",CPPFLAGS="-DEIGEN_USE_BLAS")
compile("hmm.cpp")
dyn.load(dynlib("hmm"))

## Simulate SDE
lambda <- 1
gamma <- 1
sigmaX <- 0.5
sigmaY <- 0.1
x0 <- 1
par.true <- c(x0=x0,lambda=lambda,gamma=gamma,logsX=log(sigmaX),logsY=log(sigmaY))
f <- function(x) lambda*x-gamma*x^3
g <- function(x) x*0 + sigmaX
Tsim <- 0.01; T <- 50
set.seed(1);
euler <- function(x0,f,g,tvec,dB=NULL){
    X <- numeric(length(tvec))
    X[1] <- x0
    dt <- diff(tvec)
    if(is.null(dB)) dB <- rnorm(length(dt),sd=sqrt(dt))
    for(i in 1:(length(tvec)-1))
        X[i+1] <- X[i] + f(X[i])*dt[i] + g(X[i])*dB[i]
    return(X)
}
tsim <- seq(0,T,Tsim)
Xsim <- euler(x0,f,g,tsim)
# Measure every 10th simulated value
iobs <- seq(1,length(tsim),10)
## Measurements
Y <- rnorm(length(iobs), mean = Xsim[iobs], sd = sigmaY)
plot(tsim,Xsim,type="l")
points(tsim[iobs],Y,col="red")
grid <- seq(-3,3,length=101)

data <- list(
    grid = grid,
    dt = diff(tsim[iobs])[1],
    yobs = findInterval(Y,grid)-1
    )
parameters <- list(
    lambda=0,gamma=0,logsX=0,logsY=0
    )

obj <- MakeADFun(data=data,parameters=parameters,DLL="hmm")
system.time(opt <- nlminb(obj$par,obj$fn,obj$gr))
(sdr <- sdreport(obj,opt$par))

dp <- par.true[-1] - opt$par
H <- solve(sdr$cov.fixed)
## Both less than qchisq(.95,df=4):
t(dp) %*% H %*% dp
2 * (obj$fn(par.true[-1]) - obj$fn(opt$par))
