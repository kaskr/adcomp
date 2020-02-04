## Utility
tape_print <- function(x, depth=0, dot=FALSE, DLL=getUserDLL()) {
    if (is.list(x)) x <- x$ptr
    control <- list(depth=as.integer(depth), dot=as.integer(dot))
    .Call("tmbad_print", x, control, PACKAGE=DLL)
}

src_transform <- function(obj, what=c("ADFun", "ADGrad", "ADHess"),
                          flags = "-O3") {
    if(.Platform$OS.type=="windows"){
        ## Overload tempfile
        tempfile <- function(...){
            ans <- base::tempfile(...)
            chartr("\\", "/", shortPathName(ans))
        }
    }
    what <- match.arg(what)
    DLL <- obj$env$DLL
    control <- list(method=2L)
    dll <- tempfile()
    dll.cpp <- paste0(dll, ".cpp")
    sink(dll.cpp)
    ptr <- get(what, obj$env)$ptr
    qw <- .Call("tmbad_print", ptr, control, PACKAGE = DLL)
    sink(NULL)
    ## Overload
    compile(dll.cpp, flags=flags)
    dyn.load(dynlib(dll))
    dllinfo <- getLoadedDLLs()[[basename(dll)]]
    forward_compiled <- getNativeSymbolInfo("forward",PACKAGE=dllinfo)$address
    reverse_compiled <- getNativeSymbolInfo("reverse",PACKAGE=dllinfo)$address
    .Call("TransformADFunObject",
          ptr,
          list(method="set_compiled",
               forward_compiled=forward_compiled,
               reverse_compiled=reverse_compiled
               ), PACKAGE=DLL)
}
