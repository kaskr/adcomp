library(TMB)

## Using the faster 'FFTW' library may be as simple as:
##   compile("fft.cpp", CPPFLAGS="-DEIGEN_FFTW_DEFAULT", PKG_LIBS="-lfftw3")
compile("fft.cpp")
dyn.load(dynlib("fft"))

## Data and parameters
n <- 100
d <- 0:(n-1)
d <- pmin(d, n-d)
C <- exp(-.1*d)
## x ~ MVNORM(0, C)
set.seed(1)
u <- rnorm(n)
x <- Re(fft(sqrt(fft(C)) * fft(u, TRUE)) / n)
data <- list(d=d, x=x)
parameters <- list(rho=1)

## Make a function object
obj <- MakeADFun(data, parameters, DLL="fft")

## Call function minimizer
fit <- nlminb(obj$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr
