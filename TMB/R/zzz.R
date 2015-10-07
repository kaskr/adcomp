## .First.lib <- function(lib, pkg) {
##   library.dynam("TMB", pkg, lib)
## }

.onLoad <- function(lib, pkg) {
  library.dynam("TMB", pkg, lib)
}

## .LastLib <- function(libpath)
## {
##   library.dynam.unload("TMB", libpath)
## }


