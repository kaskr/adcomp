## Copyright (C) 2013-2015 Kasper Kristensen
## License: GPL-2

##' Get or set internal configuration variables of user's DLL.
##'
##' A model compiled with the \code{TMB} C++ library has several
##' configuration variables set by default. The variables can be read
##' and modified using this function. The meaning of the variables can
##' be found in the Doxygen documentation.
##'
##' @title Get or set internal configuration variables
##' @param ... Variables to set
##' @param DLL Name of user's DLL. Auto-detected if missing.
##' @return List with current configuration
##' @examples
##' \dontrun{
##' ## Load library
##' dyn.load(dynlib("mymodel"))
##' ## Read the current settings
##' config(DLL="mymodel")
##' ## Reduce memory peak of a parallel model by creating tapes in serial
##' config(tape.parallel=0, DLL="mymodel")
##' obj <- MakeADFun(..., DLL="mymodel")
##' }
config <- function(...,DLL=getUserDLL()){
  new <- list(...)
  ## Get
  e <- new.env()
  .Call("TMBconfig",e,as.integer(1),PACKAGE=DLL)
  conf <- eapply(e,as.integer)
  ## Set
  conf[names(new)] <- new
  conf <- lapply(conf,as.integer)
  e <- local(environment(),conf)
  .Call("TMBconfig",e,as.integer(2),PACKAGE=DLL)
  ## Get
  e <- new.env()
  .Call("TMBconfig",e,as.integer(1),PACKAGE=DLL)
  conf <- eapply(e,as.integer)
  conf
}
