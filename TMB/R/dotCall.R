## Fixed R-API to .Call within MakeADFun
getParameterOrder <- function(data, parameters, ..., DLL) {
    .Call("getParameterOrder", data, parameters, new.env(), NULL, PACKAGE=DLL)
}

MakeDoubleFunObject <- function(data, parameters, reportenv, ..., DLL) {
    ans <- .Call("MakeDoubleFunObject", data, parameters, reportenv, NULL, PACKAGE=DLL)
    ans <- registerFinalizer(ans, DLL)
    ans
}
