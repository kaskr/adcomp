args <- commandArgs(TRUE)
library(TMB)
stopifnot(!dir.exists(args[2]))
dir.create(args[2], recursive=TRUE)
out <- normalizePath(args[2])
file.copy(args[1], file.path(out, 'config_test.cpp'))
setwd(out)
stopifnot(compile('config_test.cpp', framework='TMBad',
                  flags='-O1 -g0 -std=c++17', openmp=FALSE) == 0)
dyn.load(dynlib('config_test'))
dll <- 'config_test'
stopifnot(config(DLL=dll)$tmbad.sparse_hessian_batch_width == 0)
values <- list()
for (w in c(0L, 2L, 4L, 8L)) {
    config(tmbad.sparse_hessian_compress=TRUE,
           tmbad.sparse_hessian_batch_width=w, DLL=dll)
    o <- MakeADFun(list(), list(mu=.1, u=rep(.01, 11)),
                   random='u', DLL=dll, silent=TRUE)
    values[[as.character(w)]] <- unlist(lapply(c(.1, .2, -.1), function(x)
        c(o$fn(x), o$gr(x))))
}
for (v in values)
    stopifnot(all(is.finite(v)), max(abs(v - values[[1]])) < 1e-9)

# Reject invalid widths when constructing the Hessian.
config(tmbad.sparse_hessian_batch_width=3L, DLL=dll)
err <- tryCatch({
    o <- MakeADFun(list(), list(mu=.1, u=rep(.01, 11)),
                   random='u', DLL=dll, silent=TRUE)
    TMB:::sparseHessianFun(o)
    NULL
}, error=conditionMessage)
stopifnot(is.character(err),
          grepl('must be 0, 2, 4, or 8', err, fixed=TRUE))
config(tmbad.sparse_hessian_batch_width=0L, DLL=dll)
cat('PASS: default disabled; four widths agree; invalid width rejected.\n')
saveRDS(values, 'values.rds')
