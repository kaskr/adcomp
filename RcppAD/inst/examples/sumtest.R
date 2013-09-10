library(RcppAD)
dyn.load("sumtest.so")
parameters=list(x=numeric(5)+0.1 )
df <- list()

obj <- MakeADFun(data=df,parameters=parameters,random=c("^x"))
r <- obj$env$random
par <- obj$env$par
h <- obj$env$spHess(par)[r,r]
rg <- range(h-obj$env$f(par,order=2))
cat("ok\n")
print(diff(rg)<1e-12)
h

