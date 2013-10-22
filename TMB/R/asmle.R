################################################################
## Reuse code from "mle" and perhaps "mle2" (bbmle) to provide
## parameter tables, summaries, anova, profiles, confint etc.
################################################################

################################################################
## This is only a temporary solution: The methods are generally to
## inefficient - e.g. profiling. So, lets only adopt methods that
## does not require fn, gr, etc
## Original idea was to overload optim() by
## 1. Keeping objective function (because needed my mle-classes)
## 2. Attach "optim"-class attribute to output.
## HOWEVER, given the inefficiency issues this idea has been dropped.
## Instead we provide methods to convert standard optim output to
## either "mle" or "mle2" classes, only to get summary tables etc.
################################################################

## Turn optim call into mle object
as.mle <- function(opt){
  require(stats4)
  asLik <- function(fn,x){
    y <- sapply(names(x),as.name)
    cl <- as.call( c(list(as.name("c")),y ) )
    ##f <- function(...)fn(eval(cl))
    f <- function(...)stop("We do not support fn() dependent methods!")
    formals(f) <- as.list(x)
    f
  }
  if(is.matrix(opt$hessian)){
    vcov <- solve(opt$hessian)
  } else {
    npar <- c(length(opt$par))
    vcov <- matrix(NA,npar,npar)
  }
  new("mle", call = quote(mle(minuslogl = fn)), coef = opt$par, fullcoef = opt$par, 
      vcov=vcov, min = opt$value, details = opt[c("counts","convergence","message")],
      minuslogl = asLik(opt$fn,opt$par), 
      nobs = NA_integer_, method = "BFGS")
}
as.mle2 <- function(x){if(require(bbmle))as(as.mle(x),"mle2") else as.mle(x)}


if(FALSE){
  optim <- function(par,fn,...){ans <- stats::optim(par,fn,...);ans$fn <- fn;class(ans) <- "optim";ans}
  ## Test example
  ## f <- function(x).5*sum(x^2)
  ## gr <- function(x)x
  ## x0 <- c(a=3,b=4,c=5)
  ## opt <- optim(x0,f,hessian=TRUE,method="BFGS")
  print.optim <- function (x, ...) {
    print(as.mle(x),...)
  }
  summary.optim <- function (object, ...) {
    summary(as.mle2(object),...)
  }
  profile.optim <- function (fitted, ...) {
    profile(as.mle(fitted),...)
  }
  anova.optim <- function (object, ...) {
    M <- lapply(list(object, ...),as.mle2)
    li <- lapply(seq(M),function(i)substitute(M[[i]],list(i=i)))
    cl <- as.call(c( list(as.name("anova")), li  ))
    eval(cl)
  }
  confint.optim <- function (object, parm, level = 0.95, ...) {
    confint(as.mle(object), parm, level = 0.95, ...)
  }
  coef.optim <- function(x,...)x$par
  vcov.optim <- function(x,...)solve(as.matrix(x$hessian))
}
