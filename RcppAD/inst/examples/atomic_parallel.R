require(RcppAD)
dyn.load(dynlib("atomic_parallel"))
n <- 2;N <- 100
obj <- MakeADFun(type=c("ADFun","Fun"),
                 data=list(),
                 parameters=list(x=array(1:n,c(n,N))),random="x"
)
set.seed(123);x <- rnorm(n*n*n)
h <- obj$env$spHess()
.Call("InfoADFunObject",environment(obj$env$spHess)$ptr)
sum(obj$env$f(order=2)->hh) ##0.09193922
range(hh-h)
