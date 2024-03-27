## Copyright (C) 2013-2015 Kasper Kristensen
## License: GPL-2

## .First.lib <- function(lib, pkg) {
##   library.dynam("TMB", pkg, lib)
## }

## https://github.com/lme4/lme4/issues/768
## https://github.com/kaskr/adcomp/issues/387
get_abi_version <- function() {
    if (utils::packageVersion("Matrix") < "1.6-2") return(numeric_version("0"))
    Matrix::Matrix.Version()[["abi"]]
}

.Matrix.abi.build.version <- get_abi_version()

checkMatrixPackageVersion <- function(warn=TRUE) {
    cur_version <- get_abi_version()
    built_version <- .Matrix.abi.build.version
    result_ok <- cur_version == built_version
    if (!result_ok) {
        warning(
            "Package version inconsistency detected.\n",
            "TMB was built with Matrix ABI version ",
            built_version,
            "\n",
            "Current Matrix ABI version is ",
            cur_version,
            "\n",
            "Please re-install 'TMB' from source using install.packages('TMB', type = 'source') ",
            "or ask CRAN for a binary version of 'TMB' matching CRAN's 'Matrix' package"
        )
    }
    return(result_ok)
}

.onLoad <- function(lib, pkg) {
    library.dynam("TMB", pkg, lib)
    checkMatrixPackageVersion(getOption("TMB.check.Matrix", TRUE))
    ## Select AD framework (CppAD or TMBad) used by TMB::compile
    tmb.ad.framework <- getOption("tmb.ad.framework", NULL)
    if (is.null(tmb.ad.framework))
        tmb.ad.framework <- Sys.getenv("TMB_AD_FRAMEWORK", "CppAD")
    options("tmb.ad.framework" = tmb.ad.framework)
}

.onUnload <- function(libpath) {
    library.dynam.unload("TMB", libpath)
}

## .LastLib <- function(libpath)
## {
##   library.dynam.unload("TMB", libpath)
## }


