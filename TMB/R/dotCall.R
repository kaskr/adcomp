## Fixed R-API to .Call within MakeADFun
getParameterOrder <- function(data, parameters, DLL) {
    .Call("getParameterOrder", data, parameters, new.env(), NULL, PACKAGE=DLL)
}

## -----------------------------------------------------------------------------
## Constructors:

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

MakeADGradObject <- function(data, parameters, reportenv, random=NULL, f=NULL, DLL) {
    control <- list( f=f )
    if (!is.null(random))
        control$random <- as.integer(random)
    ans <- .Call("MakeADFunObject",
                 data, parameters, reportenv, control, PACKAGE=DLL)
    ans <- registerFinalizer(ans, DLL)
    ans
}

## gf   (optional) = already calculated gradient object.
## skip (optional) = index vector of parameters to skip.
MakeADHessObject <- function(data, parameters, reportenv, gf=NULL, skip=integer(0), DLL) {
    control <- list(gf=gf, skip=as.integer(skip))
    ans <- .Call("MakeADHessObject2",
                 data, parameters, reportenv, control, PACKAGE=DLL)
    ans <- registerFinalizer(ans, DLL)
    ans
}

## -----------------------------------------------------------------------------
## Evaluators

EvalDoubleFunObject <- function(Fun, theta, do_simulate=FALSE, get_reportdims=FALSE) {
    theta <- as.double(theta)
    control = list(do_simulate    = as.integer(do_simulate),
                   get_reportdims = as.integer(get_reportdims) )
    .Call("EvalDoubleFunObject", Fun$ptr, theta, control, PACKAGE=Fun$DLL)
}
