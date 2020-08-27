## Demonstrate how to pass objects back to R using the REPORT macro
library(Matrix)
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

## REPORT()
## Check objects are unchanged when passed back to R
obj <- MakeADFun(data=data, parameters=parameters, DLL="report")
(rep <- obj$report())
stopifnot( identical(rep$a, a) )
stopifnot( identical(rep$b, b) )
stopifnot( identical(rep$c, c) )
stopifnot( identical(rep$d, d) )
stopifnot( identical(rep$p, p) )
stopifnot( identical(rep$voa, list(c, c)) )
stopifnot( identical(rep$aoa, array(list(c, c), c(1,1,2,1))) )

## ADREPORT()
## Check objects are unchanged when passed back to R
sdr <- sdreport(obj, hessian.fixed=diag(1))
adrep <- as.list(sdr, "Estimate", report=TRUE)
stopifnot( all( adrep$a == a ) )
stopifnot( identical(adrep$b, b) ) ## Dimension preserved
stopifnot( identical(adrep$c, c) ) ## Dimension preserved
stopifnot( all( adrep$p == p ) )
