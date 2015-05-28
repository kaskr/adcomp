## Demonstrate how to pass objects back to R using the REPORT macro
library(TMB)

## Compile and load example
compile("report.cpp")
dyn.load(dynlib("report"))

## Generate random objects
set.seed(123)
a <- rnorm(10)
b <- array(rnorm(2*3*4),c(2,3,4))
c <- matrix(rnorm(6),2)
d <- spMatrix(10,15,i=c(1L,6L),j=c(2L,4L),x=3:4)
p <- rnorm(1)
data <- list(a=a,b=b,c=c,d=d)
parameters <- list(p=p)

## Check objects are unchanged when passed back to R
obj <- MakeADFun(data=data,parameters=parameters,DLL="report")
(rep <- obj$report())
identical(rep$a, a)
identical(rep$b, b)
identical(rep$c, c)
identical(rep$d, d)
identical(rep$p, p)
identical(rep$voa, list(c, c))
