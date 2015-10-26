## Copyright (C) 2013-2015 Kasper Kristensen
## License: GPL-2

##' Benchmark parallel templates
##'
##' By default this function will perform timings of the most critical
##' parts of an AD model, specifically
##' \enumerate{
##' \item Objective function of evaluated template.
##' \item Gradient of evaluated template.
##' \item Sparse hessian of evaluated template.
##' \item Cholesky factorization of sparse hessian.
##' }
##' (for pure fixed effect models only the first two).
##' Expressions to time can be overwritten by the user (\code{expr}).
##' A \code{plot} method is available for Parallel benchmarks.
##'
##' @title Benchmark parallel templates
##' @param obj Object from \code{MakeADFun}
##' @param n Number of replicates to obtain reliable results.
##' @param expr Optional expression to benchmark instead of default.
##' @param cores Optional vector of cores.
##' @examples
##' \dontrun{
##' runExample("linreg_parallel",thisR=TRUE)  ## Create obj
##' ben <- benchmark(obj,n=100,cores=1:4)
##' plot(ben)
##' ben <- benchmark(obj,n=10,cores=1:4,expr=expression(do.call("optim",obj)))
##' plot(ben)
##' }
benchmark <- function(obj,n=10,expr=NULL,cores=NULL){
  if(!is.null(cores)){
    return(parallelBenchmark(obj,n=n,cores=cores,expr=expr))
  }
  if(is.null(expr)){
    ## Expressions to time
    expr <- expression(
        template.likelihood = obj$env$f(order=0),
        template.gradient = obj$env$f(order=1),
        template.sparse.hessian = obj$env$spHess(random=TRUE),
        cholesky=updateCholesky(L,h)
        )
  }
  else if(is.null(names(expr)))
    names(expr) <- vapply(expr, function(.) deparse(.)[[1L]], "")
  addLoopToExpression <- function(y) substitute(for (i in seq_len(n)) { EE }, list(EE=y))
  expr <- lapply(expr, addLoopToExpression)
  if(!is.null(obj$env$random)){
    h <- obj$env$spHess() ## Avoid timing the construction
    h@x[] <- 0
    diag(h) <- 1
    L <- Cholesky(h,super=TRUE)
  } else {
    expr$template.sparse.hessian <- NULL
    expr$cholesky <- NULL
  }
  ans <- lapply(expr,function(x)system.time(eval(x)))
  ans <- do.call("rbind",ans)
  as.data.frame(ans)[c(3)]
}

## Internal helper function
parallelBenchmark <- function(obj,n,cores=1:4,...){
  ans <- lapply(cores,function(nc){
    openmp(nc)
    obj$env$retape()
    benchmark(obj,n=n,cores=NULL,...)
  })
  ans <- t(do.call("cbind",ans))
  rownames(ans) <- cores
  names(dimnames(ans)) <- c("ncores","")
  ans <- t(ans)
  class(ans) <- "parallelBenchmark"
  ans
}

##' Plot result of parallel benchmark
##'
##' @title Plot result of benchmark
##' @param x Object to plot
##' @param type Plot type
##' @param ... Further plot arguments
##' @param show Plot relative speedup or relative time?
##' @param legendpos Position of legend
##' @return NULL
##' @rdname benchmark
##' @method plot parallelBenchmark
##' @S3method plot parallelBenchmark
plot.parallelBenchmark <- function(x,type="b",...,show=c("speedup","time"),legendpos="topleft"){
  show <- match.arg(show)
  ncores <- as.numeric(colnames(x))
  x <- x[,order(ncores),drop=FALSE]
  ncores <- sort(ncores)
  if(show=="time"){
    matplot(ncores,t(x/x[,1]),ylab="Time (relative to one core)",type=type,...)
    plot(function(x)1/x,min(ncores),max(ncores),add=TRUE,col="grey")
  }
  if(show=="speedup"){
    matplot(ncores,t(x[,1]/x),type=type,ylab="Speedup (relative to one core)",...)
    abline(0,1,col="grey")
  }
  if(!is.null(legendpos)){
    n <- nrow(x)
    if(is.null(rownames(x)))rownames(x) <- 1:n
    legend(legendpos,legend=rownames(x),col=1:n,lty=1:n)
  }
  NULL
}
