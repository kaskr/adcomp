## Global options specific for RcppAD
.optionEnv <- new.env()

## ========================== Option defaults
setDefaults <- function(){
  .optionEnv$newton <- newtonDefaults()
  .optionEnv$flags <- flagsDefaults()
  .optionEnv$debug <- FALSE
}
newtonDefaults <- function(){
  list(
       trace=1,
       smartsearch=TRUE,
       mgcmax=1e60,
       tol=1e-8,
       maxit=100,
       upper = 2,
       stol=sqrt(.Machine$double.xmin),
       linesearch=FALSE,
       alpha=1,
       robustHess=FALSE,
       super=TRUE,
       silent=TRUE,
       ustep = 1, ## Start out optimistic: Newton step
       power=.5, ## decrease=function(u)const*u^power
       u0=1e-4,  ## Increase u=0 to this value  
       tol10=1e-3 ## Try to exit if last 10 iterations not improved much
       )
}
flagsDefaults <- function(){
  function(file){
    if(RcppADisParallel(file)) "-O2 -fopenmp" else NULL
  }
}

## =========================== Options extract + modify
rcppadOption <- function(x).optionEnv[[x]]
## rcppadOption <- function().GlobalEnv$.rcppadOption
newtonOption <- function(...){
  x <- list(...)
  y <- .optionEnv$newton
  ##y <- .GlobalEnv$.rcppadOption$newton
  if(length(x)==0)return(y)
  if(length(x)==1 & is.character(x[[1]]))return(y[[x[[1]]]])
  y[names(x)] <- x
  ##.GlobalEnv$.rcppadOption$newton <<- y
  .optionEnv$newton <- y
  NULL
}

