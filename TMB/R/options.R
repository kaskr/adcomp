## Global options specific for TMB
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
    ## FIXME: TMBisParallel()  does not exist ???
    if(TMBisParallel(file)) "-O2 -fopenmp" else NULL
  }
}

## =========================== Options extract + modify
tmbOption <- function(x).optionEnv[[x]]
## tmbOption <- function().GlobalEnv$.tmbOption
newtonOption <- function(...){
  x <- list(...)
  y <- .optionEnv$newton
  ##y <- .GlobalEnv$.tmbOption$newton
  if(length(x)==0)return(y)
  if(length(x)==1 & is.character(x[[1]]))return(y[[x[[1]]]])
  y[names(x)] <- x
  ##.GlobalEnv$.tmbOption$newton <<- y
  .optionEnv$newton <- y
  NULL
}

