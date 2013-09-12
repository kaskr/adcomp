config <- function(...,DLL=getUserDLL()){
  new <- list(...)
  ## Get
  e <- new.env()
  .Call("RcppADconfig",e,as.integer(1),PACKAGE=DLL)
  conf <- eapply(e,as.integer)
  ## Set
  conf[names(new)] <- new
  conf <- lapply(conf,as.integer)
  e <- local(environment(),conf)
  .Call("RcppADconfig",e,as.integer(2),PACKAGE=DLL)
  ## Get
  e <- new.env()
  .Call("RcppADconfig",e,as.integer(1),PACKAGE=DLL)
  conf <- eapply(e,as.integer)
  conf
}
