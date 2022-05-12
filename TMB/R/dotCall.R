## Fixed R-API to .Call within MakeADFun
getParameterOrder <- function(data, parameters, ..., DLL) {
    .Call("getParameterOrder", data, parameters, new.env(), NULL, PACKAGE=DLL)
}
