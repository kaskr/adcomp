## .First.lib <- function(lib, pkg) {
##   library.dynam("TMB", pkg, lib)
##   setDefaults()
## }

.onLoad <- function(lib, pkg) {
  library.dynam("TMB", pkg, lib)
  setDefaults()
}


## .LastLib <- function(libpath)
## {
##   library.dynam.unload("TMB", libpath)
## }


