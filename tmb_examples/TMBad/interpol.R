library(TMB)

## Compile and load the model
compile("interpol.cpp", framework="TMBad")
dyn.load(dynlib("interpol"))

## Smooth representation of 'volcano' data
data(volcano)
data <- list(A = volcano,
             x_range=c(0,1),
             y_range=c(0,1),
             R=2)
parameters <- list(x=0, y=0)

## Test different values of interpolation radius (R)
myf <- function(x, y, R) {
    data$R <- R
    parameters <- list(x=x, y=y)
    obj <- MakeADFun(data, parameters, DLL="interpol", silent=TRUE, type="double")
    obj$report(c(x,y))$f
}
layout(matrix(1:4,2))
x <- seq(.2,.4,length=100)
y <- seq(.2,.4,length=100)
for (R in 1:4) {
    z <- outer(x,y,myf, R=R)
    image(x,y,z)
    contour(x,y,z,add=TRUE)
    title(paste0("R=",R))
}

## Check higher order derivatives
obj <- MakeADFun(data, parameters, DLL="interpol", silent=TRUE)
library(numDeriv)
p <- c(.5, .5)
obj$fn(p)
obj$gr(p)
obj$gr(p) - grad(obj$fn, p)
obj$he(p) - jacobian(obj$gr, p)

## Optimize
opt <- optim(c(0.2,0.2), obj$fn, obj$gr, control=list(fnscale=-1))
