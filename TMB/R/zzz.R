## Copyright (C) 2013-2015 Kasper Kristensen
## License: GPL-2

## .First.lib <- function(lib, pkg) {
##   library.dynam("TMB", pkg, lib)
## }

checkMatrixPackageVersion <- function() {
    ## It is unsafe to use the TMB package with versions of 'Matrix'
    ## other than the one TMB was originally built with.
    file <- paste0(system.file(package="TMB"),"/Matrix-version")
    cur.Matrix.version <- as.character(packageVersion("Matrix"))
    if(!file.exists(file)) {
        writeLines(cur.Matrix.version, con = file)
    }
    TMB.Matrix.version <- readLines(file)
    if(!identical(TMB.Matrix.version, cur.Matrix.version)) {
        warning(
            "Package version inconsistency detected.\n",
            "TMB was built with Matrix version ",
            TMB.Matrix.version,
            "\n",
            "Current Matrix version is ",
            cur.Matrix.version,
            "\n",
            "Please re-install 'TMB' from source using install.packages('TMB', type = 'source') ",
            "or ask CRAN for a binary version of 'TMB' matching CRAN's 'Matrix' package"
        )
    }
}

.onLoad <- function(lib, pkg) {
    library.dynam("TMB", pkg, lib)
    checkMatrixPackageVersion()
}

## .LastLib <- function(libpath)
## {
##   library.dynam.unload("TMB", libpath)
## }


