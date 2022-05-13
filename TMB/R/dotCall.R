## -----------------------------------------------------------------------------
## Fixed R-API to .Call within MakeADFun
## -----------------------------------------------------------------------------

## General notes:
## - Some TMB functionality (DATA_UPDATE) implicitly assumes that 'env'
##   can be found as the enclosing environment (parent.env) of
##   'reportenv' (!). It follows that reportenv must always be passed
##   by the caller.

getParameterOrder <- function(data, parameters, reportenv, DLL) {
    control <- NULL
    .Call("getParameterOrder", data, parameters, reportenv, control, PACKAGE=DLL)
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
    ans <- .Call("MakeADGradObject",
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

EvalADFunObject <- function(ADFun, theta,
                            order=0,
                            hessiancols=NULL,
                            hessianrows=NULL,
                            sparsitypattern=FALSE,
                            rangecomponent=1,
                            rangeweight=NULL,
                            dumpstack=FALSE,
                            doforward=TRUE,
                            set_tail=FALSE,
                            keepx=NULL,
                            keepy=NULL,
                            data_changed=FALSE) {
    if (!is.null(rangeweight))
        rangeweight <- as.double(rangeweight)
    control <- list(order=as.integer(order),
                    hessiancols=as.integer(hessiancols),
                    hessianrows=as.integer(hessianrows),
                    sparsitypattern=as.integer(sparsitypattern),
                    rangecomponent=as.integer(rangecomponent),
                    rangeweight=rangeweight,
                    dumpstack=as.integer(dumpstack),
                    doforward=as.integer(doforward),
                    set_tail = as.integer(set_tail),
                    keepx=as.integer(keepx),
                    keepy=as.integer(keepy),
                    data_changed = as.integer(data_changed) )
    .Call("EvalADFunObject", ADFun$ptr, theta, control, PACKAGE=ADFun$DLL)
}
