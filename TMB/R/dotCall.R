## Fixed R-API to .Call within MakeADFun
getParameterOrder <- function(data, parameters, DLL) {
    .Call("getParameterOrder", data, parameters, new.env(), NULL, PACKAGE=DLL)
}

## Constructor input:
## - data
## - parameters
## - environment
## - control

MakeDoubleFunObject <- function(data, parameters, reportenv, DLL) {
    control <- NULL
    ans <- .Call("MakeDoubleFunObject",
                 data, parameters, reportenv, control, PACKAGE=DLL)
    ans <- registerFinalizer(ans, DLL)
    ans
}

MakeADFunObject <- function(data, parameters, reportenv, ADreport=FALSE, DLL) {
    control <- list( report = as.integer(ADreport) )
    ans <- .Call("MakeADFunObject",
                 data, parameters, reportenv, control, PACKAGE=DLL)
    ans <- registerFinalizer(ans, DLL)
    ans
}
