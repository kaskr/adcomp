TransformADFunObject <- function(ADFun,
                                 method,
                                 ...
                                 ) {
    method <- as.character(method)
    ans <- .Call("TransformADFunObject",
                 f = ADFun$ptr,
                 control = list(method = method, ...),
                 PACKAGE = ADFun$DLL)
    if (method == "copy") {
        ans <- registerFinalizer(ans, ADFun$DLL)
    }
    ans
}

## Not really a transform but didn't know where else to put it...
readNodeInputs <- function(ADFun, node) {
    node <- as.double(node)
    TransformADFunObject(ADFun, "readNodeInputs", node=node)
}

## Utility
tape_print <- function(x, depth=0, method="tape", DLL=getUserDLL(), ...) {
    if (is.list(x)) x <- x$ptr
    control <- list(depth=as.integer(depth), method=as.character(method), ...)
    .Call("tmbad_print", x, control, PACKAGE=DLL)
}

op_table <- function(ADFun, name=TRUE, address=FALSE, input_size=FALSE, output_size=FALSE) {
    ntapes <- tape_print(ADFun, method="num_tapes", DLL=ADFun$DLL, i=as.integer(0))
    ntapes <- max(1, ntapes)
    f <- function(i)tape_print(ADFun$ptr, method="op", DLL=ADFun$DLL, i=as.integer(i),
                                     name=as.integer(name),
                                     address=as.integer(address),
                                     input_size=as.integer(input_size),
                                     output_size=as.integer(output_size))
    g <- function(i)data.frame(tape=i, opname=f(i), stringsAsFactors=FALSE)
    df <- do.call("rbind", lapply(seq_len(ntapes) - 1L, g))
    table(opname = df$opname, tape = df$tape)
}

src_transform <- function(ADFun,
                          flags = "-O3", ..., perm=TRUE) {
    if(.Platform$OS.type=="windows"){
        ## Overload tempfile
        tempfile <- function(...){
            ans <- base::tempfile(...)
            chartr("\\", "/", shortPathName(ans))
        }
    }
    ntapes <- tape_print(ADFun, method="num_tapes",
                               DLL=ADFun$DLL,
                               i=as.integer(0))
    ntapes <- max(1, ntapes)
    tapes <- seq.int(from=0, length.out=ntapes)
    control <- list(method="src")
    dll <- tempfile(fileext=paste0("_",tapes))
    dll.cpp <- paste0(dll, ".cpp")
    ## Reorder graph
    if (perm) {
        TransformADFunObject(
                  ADFun,
                  method="reorder_sub_expressions",
                  random_order=integer(0),
                  max_period_size=1024L)
    }
    ## Write redefs
    forward <- paste0("forward", tapes)
    reverse <- paste0("reverse", tapes)
    redef <- function(i) {
        cat("#define forward", forward[i+1], "\n")
        cat("#define reverse", reverse[i+1], "\n")
    }
    ## Write source code
    for (i in tapes) {
        control$i <- i
        sink(dll.cpp[i+1]); redef(i); out <- .Call("tmbad_print", ADFun$ptr, control, PACKAGE = ADFun$DLL); sink(NULL)
    }
    ## Overload
    compile(dll.cpp, flags=flags, ..., libtmb=FALSE)
    dyn.load(dynlib(dll)[1])
    dllinfo <- getLoadedDLLs()[[basename(dll[1])]]
    forward_compiled <-
        lapply(forward, function(x)getNativeSymbolInfo(x,PACKAGE=dllinfo)$address)
    reverse_compiled <-
        lapply(reverse, function(x)getNativeSymbolInfo(x,PACKAGE=dllinfo)$address)
    TransformADFunObject(
              ADFun,
              method="set_compiled",
              forward_compiled=forward_compiled,
              reverse_compiled=reverse_compiled)
    ## Unload compiled code when no longer needed
    finalizer <- function(ptr) {
        dyn.unload(dynlib(dll[1]))
        file.remove(dynlib(dll[1]))
        file.remove(paste0(dll, ".o"))
        file.remove(dll.cpp)
    }
    reg.finalizer(ADFun$ptr, finalizer)
    NULL
}
