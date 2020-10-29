TransformADFunObject <- function(ADFun,
                                 method,
                                 ...
                                 ) {
    .Call("TransformADFunObject",
          f = ADFun$ptr,
          control = list(method = as.character(method), ...),
          PACKAGE = ADFun$DLL)
}

## Utility
tape_print <- function(x, depth=0, method="tape", DLL=getUserDLL(), ...) {
    if (is.list(x)) x <- x$ptr
    control <- list(depth=as.integer(depth), method=as.character(method), ...)
    .Call("tmbad_print", x, control, PACKAGE=DLL)
}

op_table <- function(ADFun) {
    ntapes <- TMB:::tape_print(ADFun, method="num_tapes", DLL=ADFun$DLL, i=as.integer(0))
    ntapes <- max(1, ntapes)
    f <- function(i)TMB:::tape_print(ADFun$ptr, method="opname", DLL=ADFun$DLL, i=as.integer(i))
    g <- function(i)data.frame(tape=i, opname=f(i), stringsAsFactors=FALSE)
    df <- do.call("rbind", lapply(seq_len(ntapes) - 1L, g))
    table(opname = df$opname, tape = df$tape)
}

src_transform <- function(obj, what=c("ADFun", "ADGrad", "ADHess"),
                          flags = "-O3", perm=TRUE) {
    if(.Platform$OS.type=="windows"){
        ## Overload tempfile
        tempfile <- function(...){
            ans <- base::tempfile(...)
            chartr("\\", "/", shortPathName(ans))
        }
    }
    what <- match.arg(what)
    DLL <- obj$env$DLL
    control <- list(method="src")
    dll <- tempfile()
    dll.cpp <- paste0(dll, ".cpp")
    ptr <- get(what, obj$env)$ptr
    ## Reorder graph
    if (perm) {
        .Call("TransformADFunObject",
              ptr,
              list(random_order=integer(0),
                   max_period_size=1024L,
                   method="reorder_sub_expressions"), PACKAGE=DLL)
    }
    ## Write source code
    sink(dll.cpp); out <- .Call("tmbad_print", ptr, control, PACKAGE = DLL); sink(NULL)
    ## Overload
    compile(dll.cpp, flags=flags, libtmb=FALSE)
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
    ## Unload compiled code when no longer needed
    finalizer <- function(ptr) {
        dyn.unload(dynlib(dll))
        file.remove(dynlib(dll))
        file.remove(paste0(dll, ".o"))
        file.remove(dll.cpp)
    }
    reg.finalizer(ptr, finalizer)
    NULL
}
