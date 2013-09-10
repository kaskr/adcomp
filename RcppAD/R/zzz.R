## .First.lib <- function(lib, pkg) {
##   library.dynam("RcppAD", pkg, lib)
##   setDefaults()
## }

.onLoad <- function(lib, pkg) {
  library.dynam("RcppAD", pkg, lib)
  setDefaults()
}


## .LastLib <- function(libpath)
## {
##   library.dynam.unload("RcppAD", libpath)
## }


