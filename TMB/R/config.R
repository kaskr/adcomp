## Copyright (C) 2013-2015 Kasper Kristensen
## License: GPL-2

##' Get or set internal configuration variables of user's DLL.
##'
##' A model compiled with the \code{TMB} C++ library has several
##' configuration variables set by default. The variables can be read
##' and modified using this function. The meaning of the variables can
##' be found in the Doxygen documentation.
##'
##' @section Sparse Hessian memory with TMBad:
##' For models compiled with \code{framework="TMBad"},
##' \code{tmbad.sparse_hessian_compress=TRUE} can reduce the memory needed
##' to represent the sparse random-effect Hessian used in the Laplace
##' approximation. Set this option after loading the model DLL and before
##' creating the objective with \code{\link{MakeADFun}}. It does not
##' change which parameters are treated as random effects.
##'
##' Compression is disabled by default. Its performance depends on the model:
##' lower memory use can come with slower objective or gradient evaluations.
##' Compare both first-use and repeated evaluations at identical parameter
##' values before choosing a setting. When comparing settings, create a new
##' objective for each setting; changing the configuration does not rebuild
##' derivative tapes that have already been created.
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
##' ## For a DLL compiled with framework="TMBad", try Hessian compression
##' config(tmbad.sparse_hessian_compress=TRUE, DLL="mymodel")
##' compressed <- MakeADFun(data, parameters, random="u", DLL="mymodel")
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
  ## Backwards compatibility (don't break if members unavailable)
  if (is.null(conf$autopar))  conf$autopar  <- FALSE
  if (is.null(conf$nthreads)) conf$nthreads <- 1
  conf
}
