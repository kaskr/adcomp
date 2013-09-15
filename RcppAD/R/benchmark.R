benchmark <- function(obj,n=10){
  h <- obj$env$spHess() ## Avoid timing the construction
  h@x[] <- 0
  diag(h) <- 1
  L <- Cholesky(h,super=TRUE)
  expr <- expression(
      template.likelihood = for(i in seq(length=n))obj$env$f(order=0),
      template.gradient = for(i in seq(length=n))obj$env$f(order=1),
      template.sparse.hessian = for(i in seq(length=n))obj$env$spHess(),
      cholesky=for(i in seq(length=n))updateCholesky(L,h)
      )
  ans <- lapply(expr,function(x)system.time(eval(x)))
  ans <- do.call("rbind",ans)
  as.data.frame(ans)[c(3)]
}

benchmarkParallel <- function(obj,n,cores=1:4){
  ans <- lapply(cores,function(nc){
    openmp(nc)
    obj$env$retape()
    benchmark(obj,n=n)
  })
  ans <- t(do.call("cbind",ans))
  rownames(ans) <- cores
  names(dimnames(ans)) <- c("ncores","")
  ans <- t(ans)
  class(ans) <- "benchmarkParallel"
  ans
}

plot.benchmarkParallel <- function(x,...){
  ncores <- as.numeric(colnames(x))
  matplot(ncores,t(x/x[,1]),ylab="Time (relative to one core)",...)
}
