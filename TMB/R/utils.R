## Get *subset* of inverse of sparse matrix Q.
## Subset is guarantied to contain pattern of Q.
## Q: Sparse positive definite matrix.
## L: Cholesky factor of Q.
## diag: Return just diagonal of inverse ?
solveSubset <- function(Q,
                        L = Matrix::Cholesky(Q, super=TRUE, perm=TRUE),
                        diag = FALSE) {
    stopifnot( is(L, "dCHMsuper") )
    invQ <- .Call("tmb_invQ", L, PACKAGE = "TMB")
    iperm <- Matrix::invPerm(L@perm + 1L)
    if (diag) {
        invQ <- diag(invQ)[iperm]
    } else {
        invQ <- invQ[iperm, iperm, drop=FALSE]
    }
    invQ
}

## Get information on ADFun object pointer
info <- function(ADFun, DLL = getUserDLL()) {
    ptr <- ADFun$ptr
    DLL <- ADFun$DLL
    ok <- is(ptr, "externalptr") && !isNullPointer(ptr)
    if (!ok) stop("'ptr' is not a valid external pointer")
    .Call("InfoADFunObject", ptr, PACKAGE=DLL)
}
