## Utility
tape_print <- function(x, depth=0, dot=FALSE, DLL=getUserDLL()) {
    if (is.list(x)) x <- x$ptr
    control <- list(depth=as.integer(depth), dot=as.integer(dot))
    .Call("tmbad_print", x, control, PACKAGE=DLL)
}
