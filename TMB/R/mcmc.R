## Copyright (C) 2015 Cole Monnahan
## License: GPL-2

#' [BETA VERSION] Draw samples from the posterior of a TMB model using a
#' specified MCMC algorithm.
#'
#' @details This function is a top-level wrapper designed specifically to
#'   work with TMB models. There are several MCMC algorithms available for
#'   use. The user is responsible for specifying the model properly
#'   (priors, starting values, desired parameters fixed, etc.), as well as
#'   assessing the convergence and validity of the resulting samples (e.g.,
#'   through the \code{coda} package) before making inference.
#' @title MCMC sampling of TMB models
#' @author Cole Monnahan
#' @param obj A TMB model object.
#' @param iter The number of (dependent) samples to draw.
#' @param init A list of lists containing the initial parameter vectors,
#'   one for each chain. It is strongly recommended to initialize multiple
#'   chains from dispersed points. The default of NULL signifies to use the
#'   starting values present in the model (i.e., \code{obj$par}), but this
#'   will apply to all chains.
#' @param covar An optional covariance matrix which can be used to improve
#'   the efficiency of sampling. The lower Cholesky decomposition of this
#'   matrix is used to transform the parameter space. If the posterior is
#'   approximately multivariate normal and \code{covar} approximates the
#'   covariance, then the transformed parameter space will be close to
#'   multivariate standard normal. In this case the algorithm will be more
#'   efficient, but there will be overhead in the matrix calculations which
#'   need to be done at each step. The default of NULL specifies to not do
#'   this transformation.
#' @param algorithm A string specifiying an algorithm. Currently supported
#'   are: \itemize{ \item{"RWM"}{the random walk Metropolis sampler}
#'   \item{"HMC"}{the Hamiltonian sampler (see Neal 2011)}
#'   \item{"NUTS"}{the No-U-Turn sampler (see Hoffman and Gelman 2014)} }
#'   These algorithms require different arguments; see their help files for
#'   more information.
#' @param thin The thinning rate to apply to samples. The default of 1
#'   means to not thin (keep all samples). Thinning usually is not
#'   necessary for NUTS but in some cases may be useful.
#' @param ... Further arguments to be passed to the algorithm. See help
#'   files for the samplers for further arguments.
#' @return A list containing the samples and properties of the sampler
#'   useful for diagnosing behavior and efficiency.
#' @example inst/examples/mcmc_examples.R
#' @seealso \code{\link{run_mcmc.hmc}}, \code{\link{run_mcmc.nuts}},
#'   \code{\link{run_mcmc.rwm}}
run_mcmc <- function(obj, iter, algorithm="NUTS", chains=1, init=NULL,
                     covar=NULL, lower=NULL, upper=NULL, thin=1, ...){
  ## Argument checking
  if(is.null(init)){
    if(chains>1) warning('Using same inits for each chain -- strongly recommended to use dispersed inits')
    init <- rep(list(obj$par), times=chains)
  } else if(length(init) != chains){
    stop("Length of init does not equal number of chains.")
  } else if(any(unlist(lapply(init, function(x) length(x) != length(obj$par))))){
    stop("Initial parameter vector is wrong length")
  }
  algorithm <- match.arg(algorithm, choices=c("NUTS", "RWM", "HMC"))
  thin <- floor(thin)
  stopifnot(thin >=1)
  stopifnot(chains >= 1)
  if(iter < 10 | !is.numeric(iter)) stop("iter must be > 10")
  obj$env$beSilent()                  # silence console output

  ## Parameter constraints, if provided, require the fn and gr functions to
  ## be modified to account for differents in volume. There are four cases:
  ## no constraints, bounded below, bounded above, or both (box
  ## constraint).
  bounded <- !(is.null(lower) & is.null(upper))
  if(bounded){
    if(is.null(lower)) lower <- rep(-Inf, len=length(upper))
    if(is.null(upper)) upper <- rep(Inf, len=length(lower))
    cases <- .transform.cases(lower, upper)
    fn <- function(y){
      x <- .transform(y, lower, upper, cases)
      scales <- .transform.grad(y, lower, upper, cases)
      -obj$fn(x) + sum(log(scales))
    }
    gr <- function(y){
      x <- .transform(y, lower, upper, cases)
      scales <- .transform.grad(y, lower, upper, cases)
      scales2 <- .transform.grad2(y, lower, upper, cases)
      -as.vector(obj$gr(x))*scales + scales2
    }
    init <- lapply(init, FUN=.transform.inv, a=lower, b=upper, cases=cases)
  } else {
    fn <- function(x) -obj$fn(x)
    gr <- function(x) -as.vector(obj$gr(x))
  }

  ## Make parameter names unique if vectors exist
  par.names <- names(obj$par)
  par.names <- as.vector((unlist(sapply(unique(par.names), function(x){
    temp <- par.names[par.names==x]
    if(length(temp)>1) paste0(temp,'[',1:length(temp),']') else temp
    }))))

  ## Select and run the chain.
  if(algorithm=="HMC"){
    mcmc.out <- lapply(1:chains, function(i)
      run_mcmc.hmc(iter=iter, fn=fn, gr=gr, init=init[[i]],
                   covar=covar, chain=i, thin=thin, ...))
  } else if(algorithm=="NUTS"){
    mcmc.out <- lapply(1:chains, function(i)
      run_mcmc.nuts(iter=iter, fn=fn, gr=gr, init=init[[i]],
                    covar=covar, chain=i, thin=thin, ...))
  } else if(algorithm=="RWM")
    mcmc.out <- lapply(1:chains, function(i)
      run_mcmc.rwm(iter=iter, fn=fn, init=init[[i]], covar=covar,
                  thin=thin, ...))

  ## Clean up returned output
  samples <-  array(NA, dim=c(nrow(mcmc.out[[1]]$par), chains, 1+length(par.names)),
                    dimnames=list(NULL, NULL, c(par.names,'lp__')))
  for(i in 1:chains){
    if(bounded){
      temp <- mcmc.out[[i]]$par
      temp[,-ncol(temp)] <-
        t(apply(temp[,-ncol(temp)], 1, function(x)
          .transform(x, lower, upper, cases)))
      samples[,i,] <- temp
    } else {
      samples[,i,] <- mcmc.out[[i]]$par
    }
  }
  sampler_params <- lapply(mcmc.out, function(x) x$sampler_params)
  time.warmup <- unlist(lapply(mcmc.out, function(x) as.numeric(x$time.warmup)))
  time.total <- unlist(lapply(mcmc.out, function(x) as.numeric(x$time.total)))
  result <- list(samples=samples, sampler_params=sampler_params,
                 time.warmup=time.warmup, time.total=time.total,
                 algorithm=algorithm, warmup=mcmc.out[[1]]$warmup,
                 model=obj$env$DLL)
  if(algorithm=="NUTS") result$max_treedepth <- mcmc.out[[1]]$max_treedepth
  return(invisible(result))
}
#' Determine which transformation is to be used for each parameter.
#'
#' @details The 4 cases are (0) none, (1) lower only [a,Inf], (2) upper
#'   only [-Inf, b], and (3) both lower and upper [a,b]. Each case requires
#'   a different set of transformation functions.
#' @param lower Vector of lower bounds, potentially -infinity for some
#' @param upper Vector of upper bounds, potentially infinity for some
#' @return Vector of cases, in 0:3, to be used in all other transformation
#'   functions. Error checking is only done here, not in other functions.
#' @seealso \code{\link{.transform}}, \code{\link{.transform.inv}},
#'   \code{\link{.transform.grad}}, \code{\link{.transform.grad2}}
#'
.transform.cases <- function(lower, upper){
  if(length(lower) != length(upper))
    stop("Lengths of lower and upper do not match")
  if(any(is.na(c(lower, upper))) | any(is.nan(c(lower, upper))))
    stop("Bounds must be finite or -Inf/Inf -- NA and NaN not allowed")
  if(any(lower >= upper))
    stop("Lower bound >= upper bound")
  cases <- rep(NA, length(lower))
  cases[!is.finite(lower) & !is.finite(upper)] <- 0
  cases[is.finite(lower) & !is.finite(upper)] <- 1
  cases[!is.finite(lower) & is.finite(upper)] <- 2
  cases[is.finite(lower) & is.finite(upper)] <- 3
  if(any(is.na(cases)))
    stop("Something unexpected went wrong determining the bounding functions.
 Check lower and upper.")
  return(cases)
}
#' This function returns the transformed variable, x=f(y).
#'
.transform <- function(y, a, b, cases){
  x <- sapply(1:length(y), function(i) {
    if(cases[i]==0) return(y[i])
    else if(cases[i]==1) return(exp(y[i])+a[i])
    else if(cases[i]==2) return(b[i]-exp(y[i]))
    else if(cases[i]==3) return(a[i]+(b[i]-a[i])/(1+exp(-y[i])))
  })
  return(x)
}
#' The inverse of the transformation, y=f-1(x).
#'
.transform.inv <- function(x, a, b, cases){
  if(any(x<a) | any(x>b)) stop("x outside limits provided -- not meaningful")
  y <- sapply(1:length(x), function(i) {
    if(cases[i]==0) return(x[i])
    else if(cases[i]==1) return(log(x[i]-a[i]))
    else if(cases[i]==2) return(log(b[i]-x[i]))
    else if(cases[i]==3) return(-log( (b[i]-x[i])/(x[i]-a[i]) ))
  })
  return(y)
}

#' The absolute value of the derivative of transformation.
#'
.transform.grad <- function(y, a, b, case){
  x <- sapply(1:length(y), function(i) {
    if(case[i]==0) return(1)
    else if(case[i]==1) return(exp(y[i]))
    else if(case[i]==2) return(exp(y[i]))
    else if(case[i]==3) return((b[i]-a[i])*exp(-y[i])/(1+exp(-y[i]))^2)
  })
  return(x)
}
#' The derivative of the log of the derivative of the transformation. I.e.,
#' d/dy[log(.transform.grad(y,a,b))].
#'
.transform.grad2 <- function(y, a, b, case){
  x <- sapply(1:length(y), function(i) {
    if(case[i]==0) return(0)
    else if(case[i]==1) return(1)
    else if(case[i]==2) return(1)
    else if(case[i]==3) return(-1+2*exp(y[i])/(1+exp(y[i])))
  })
  return(x)
}

#' [BETA VERSION] Draw MCMC samples from a model posterior using a
#' Random Walk Metropolis (RWM) sampler.
#'
#' @param iter The number of samples to return.
#' @param fn A function that returns the log of the posterior density.
#' @param init A vector of initial parameter values.
#' @param diagnostic Whether to return a list of diagnostic metrics about
#' the chain. Useful for assessing efficiency and tuning chain.
#' @details This algorithm does not yet contain adaptation of \code{alpha}
#' so some trial and error may be required for efficient sampling.
#' @param covar An optional covariance matrix which can be used to improve
#' the efficiency of sampling. The lower Cholesky decomposition of this
#' matrix is used to transform the parameter space. If the posterior is
#' approximately multivariate normal and \code{covar} approximates the
#' covariance, then the transformed parameter space will be close to
#' multivariate standard normal. In this case the algorithm will be more
#' efficient, but there will be overhead in the matrix calculations which
#' need to be done at each step. The default of NULL specifies to not do
#' this transformation.
#' @param alpha The amount to scale the proposal, i.e,
#' Xnew=Xcur+alpha*Xproposed where Xproposed is generated from a mean-zero
#' multivariate normal. Varying \code{alpha} varies the acceptance rate.
#' @return If \code{diagnostic} is FALSE (default), returns a matrix of
#' \code{iter} samples from the posterior. Otherwise returns a list
#' containing samples ('par'), proposed samples ('par.proposed'), vector of
#' which proposals were accepted ('accepted'), and the total function calls
#' ('n.calls'), which for this algorithm is \code{iter}
#' @seealso \code{\link{run_mcmc}}, \code{\link{run_mcmc.nuts}}, \code{\link{run_mcmc.hmc}}
run_mcmc.rwm <- function(iter, fn, init, alpha=1, chain=1,
                         warmup=floor(iter/2), covar=NULL, thin=1){
  lp <- accepted <- rep(0, length=iter)
  n.params <- length(init)
  theta.out <- matrix(NA, nrow=iter, ncol=n.params)
  ## If using covariance matrix and Cholesky decomposition, redefine
  ## these functions to include this transformation. The algorithm will
  ## work in the transformed space.
  if(!is.null(covar)){
    fn2 <- function(theta) fn(chd %*% theta)
    chd <- t(chol(covar))               # lower triangular Cholesky decomp.
    chd.inv <- solve(chd)               # inverse
    theta.cur <- chd.inv %*% init
  } else {
    fn2 <- fn
    theta.cur <- init
  }
  fn.cur <- fn2(theta.cur)
  time.start <- Sys.time()
  message('')
  message(paste('Starting RWM at', time.start))
  for(m in 1:iter){
    ## generate proposal
    theta.new <- theta.cur + alpha*rnorm(n=n.params, mean=0, sd=1)
    fn.new <- fn2(theta.new)
    if(log(runif(1))< fn.new-fn.cur){
      ## accept
      accepted[m] <- 1
      theta.cur <- theta.out[m,] <- theta.new
      fn.cur <- fn.new
    } else {
      ## do not accept
      theta.out[m,] <- theta.cur
    }
    lp[m] <- fn.cur
    if(m==warmup) time.warmup <- difftime(Sys.time(), time.start, units='secs')
    .print.mcmc.progress(m, iter, warmup, chain)
  } ## end of MCMC loop

  ## Back transform parameters if covar is used
  if(!is.null(covar)) {
    theta.out <- t(apply(theta.out, 1, function(x) chd %*% x))
  }
  theta.out <- cbind(theta.out, lp)
  theta.out <- theta.out[seq(1, nrow(theta.out), by=thin),]
  message(paste0("Final acceptance ratio=", sprintf("%.2f", mean(accepted[-(1:warmup)]))))
  time.total <- difftime(Sys.time(), time.start, units='secs')
  .print.mcmc.timing(time.warmup=time.warmup, time.total=time.total)
  return(list(par=theta.out, sampler_params=NULL,
              time.total=time.total, time.warmup=time.warmup, warmup=warmup/thin))
}

#' [BETA VERSION] Draw MCMC samples from a model posterior using a static
#' HMC sampler.
#' @details This function implements algorithm 5 of Hoffman and Gelman
#'   (2014), which includes adaptive step sizes (\code{eps}) via an
#'   algorithm called dual averaging.
#' @param iter The number of samples to return.
#' @param L The number of leapfrog steps to take. The NUTS algorithm does
#'   not require this as an input. If \code{L=1} this function will perform
#'   Langevin sampling. In some contexts \code{L} can roughly be thought of
#'   as a thinning rate.
#' @param eps The step size. If a numeric value is passed, it will be used
#'   throughout the entire chain. A \code{NULL} value will initiate
#'   sampler_params of \code{eps} using the dual averaging algorithm during
#'   the first \code{warmup} steps.
#' @param warmup How many iterations to use for a warmup, in which the step
#'   size will be adapted. The default is \code{warmup=iter/2}.
#' @param adapt_delta The target acceptance rate if using apative
#'   \code{eps}. Defaults to 0.8.
#' @param fn A function that returns the log of the posterior density.
#' @param gr A function that returns a vector of gradients of the log of
#'   the posterior density (same as \code{fn}).
#' @param init A vector of initial parameter values.
#' @param covar An optional covariance (mass) matrix which can be used to
#'   improve the efficiency of sampling. The lower Cholesky decomposition
#'   of this matrix is used to transform the parameter space. If the
#'   posterior is approximately multivariate normal and \code{covar}
#'   approximates the covariance, then the transformed parameter space will
#'   be close to multivariate standard normal. In this case the algorithm
#'   will be more efficient, but there will be overhead in the matrix
#'   calculations which need to be done at each step. The default of NULL
#'   specifies to not do this transformation and use a unit diagonal
#'   matrix.
#' @param chain The MCMC chain to run. Only used for bookkeeping at the
#'   moment.
#' @references \itemize{ \item{Neal, R. M. (2011). MCMC using Hamiltonian
#'   dynamics. Handbook of Markov Chain Monte Carlo.}  \item{Hoffman and
#'   Gelman (2014). The No-U-Turn sampler: Adaptively setting path lengths
#'   in Hamiltonian Monte Carlo. J. Mach. Learn. Res.  15:1593-1623.}  }
#' @seealso \code{\link{run_mcmc}}, \code{\link{run_mcmc.nuts}},
#'   \code{\link{run_mcmc.rwm}}
#' @return A list containing samples ('par') and algorithm details such as
#'   step size adaptation and acceptance probabilities per iteration
#'   ('sampler_params').
run_mcmc.hmc <- function(iter, fn, gr, init, L, eps=NULL, covar=NULL,
                         adapt_delta=0.8, warmup=floor(iter/2),
                         chain=1,thin=1){
  warning("NUTS should be prefered to sHMC except in rare, specific cases")
  ## If using covariance matrix and Cholesky decomposition, redefine
  ## these functions to include this transformation. The algorithm will
  ## work in the transformed space
  if(!is.null(covar)){
    fn2 <- function(theta) fn(chd %*% theta)
    gr2 <- function(theta) as.vector( t( gr(chd %*% theta) ) %*% chd )
    chd <- t(chol(covar))               # lower triangular Cholesky decomp.
    chd.inv <- solve(chd)               # inverse
    theta.cur <- chd.inv %*% init
  } else {
    fn2 <- fn; gr2 <- gr
    theta.cur <- init
  }
  accepted <- divergence <- lp <- rep(NA, iter)
  theta.out <- matrix(NA, nrow=iter, ncol=length(init))
  sampler_params <- matrix(numeric(0), nrow=iter, ncol=4, # holds DA info by iteration
                       dimnames=list(NULL, c("accept_stat__",
                                                "stepsize__", "int_time__", "energy__")))
  ## A NULL value for eps signifies to use the dual averaging algorithm
  useDA <- is.null(eps)
  if(useDA){
    ## Initialize the dual-averaging algorithm.
    epsvec <- Hbar <- epsbar <- rep(NA, length=warmup+1)
    eps <- epsvec[1] <- epsbar[1] <-
      .find.epsilon(theta=theta.cur, fn=fn2, gr=gr2, eps=.1, verbose=FALSE)
    mu <- log(10*eps)
    Hbar[1] <- 0; gamma <- 0.05; t0 <- 10; kappa <- 0.75
  } else {
    ## dummy values to return
    epsvec <- epsbar <- Hbar <- NULL
  }
  ## Start of MCMC chain
  time.start <- Sys.time()
  message('')
  message(paste('Starting static HMC at', time.start))
  for(m in 1:iter){
    ## Jitter step size to mitigate potential negative autocorrelations,
    ## only once fixed though
    if(useDA & m > warmup) eps <- eps*runif(1,.9,1.1)
    r.cur <- r.new <- rnorm(length(init),0,1)
    theta.new <- theta.cur
    theta.leapfrog <- matrix(NA, nrow=L, ncol=length(theta.cur))
    r.leapfrog <- matrix(NA, nrow=L, ncol=length(theta.cur))
    ## Make a half step for first iteration
    r.new <- r.new+eps*gr2(theta.new)/2
    for(i in 1:L){
      theta.leapfrog[i,] <- theta.new
      r.leapfrog[i,] <- r.new
      theta.new <- theta.new+eps*r.new
      ## Full step except at end
      if(i!=L) r.new <- r.new+eps*gr2(theta.new)
      ## If divergence, stop trajectory earlier to save computation
      if(any(!is.finite(r.new)) | any(!is.finite(theta.new))) break
    }
    ## half step for momentum at the end
    r.new <- r.new+eps*gr2(theta.new)/2
    logalpha <- -fn2(theta.cur)+fn2(theta.new)+ sum(r.cur^2)/2-sum(r.new^2)/2
    ## Numerical divergence is registered as a NaN above. In this case we
    ## want to reject the proposal, mark the divergence, and adjust the
    ## step size down if still adapting (see below).
    if(!is.finite(logalpha)){
      divergence[m] <- 1
      logalpha <- -Inf
    } else {
      divergence[m] <- 0
    }
    ## Test whether to accept or reject proposed state
    if(log(runif(1)) < logalpha){
      ## accept the proposed state
      theta.cur <- theta.new
      accepted[m] <- TRUE
    } else {
      ## otherwise reject it and stay there
      accepted[m] <- FALSE
    }
    theta.out[m,] <- theta.cur
    lp[m] <- fn(theta.cur)
    if(useDA){
      ## Do the adapting of eps.
      if(m <= warmup){
        Hbar[m+1] <-
          (1-1/(m+t0))*Hbar[m] + (adapt_delta-min(1,exp(logalpha)))/(m+t0)
        logeps <- mu-sqrt(m)*Hbar[m+1]/gamma
        epsvec[m+1] <- exp(logeps)
        logepsbar <- m^(-kappa)*logeps + (1-m^(-kappa))*log(epsbar[m])
        epsbar[m+1] <- exp(logepsbar)
        eps <- epsvec[m+1]
      } else {
        eps <- epsbar[warmup]
      }
    }
    ## Save adaptation info.
    sampler_params[m,] <- c(min(1,exp(logalpha)), eps, eps*L, fn2(theta.cur))
    if(m==warmup) time.warmup <- difftime(Sys.time(), time.start, units='secs')
    .print.mcmc.progress(m, iter, warmup, chain)
  } ## end of MCMC loop
  ## Back transform parameters if covar is used
  if(!is.null(covar)) {
    theta.out <- t(apply(theta.out, 1, function(x) chd %*% x))
  }
  theta.out <- cbind(theta.out, lp)
  theta.out <- theta.out[seq(1, nrow(theta.out), by=thin),]
  sampler_params <- sampler_params[seq(1, nrow(sampler_params), by=thin),]
  if(sum(divergence[-(1:warmup)])>0)
    message(paste0("There were ", sum(divergence[-(1:warmup)]),
                   " divergent transitions after warmup"))
  message(paste0("Final acceptance ratio=", sprintf("%.2f", mean(accepted[-(1:warmup)])),
                 " and target is ", adapt_delta))
  if(useDA) message(paste0("Final step size=", round(epsbar[warmup], 3),
                           "; after ", warmup, " warmup iterations"))
  time.total <- difftime(Sys.time(), time.start, units='secs')
  .print.mcmc.timing(time.warmup=time.warmup, time.total=time.total)
  return(list(par=theta.out, sampler_params=sampler_params,
              time.total=time.total, time.warmup=time.warmup, warmup=warmup/thin))
}


#' [BETA VERSION] Draw MCMC samples from a model posterior using the
#' No-U-Turn (NUTS) sampler with dual averaging.
#'
#' @details This function implements algorithm 6 of Hoffman and Gelman
#'   (2014), which includes adaptive step sizes (\code{eps}) via an
#'   algorithm called dual averaging. In theory neither the step length nor
#'   step size needs to be input by the user to obtain efficient sampling
#'   from the posterior.
#' @param iter The number of samples to return.
#' @param eps The length of the leapfrog steps. If a numeric value is
#'   passed, it will be used throughout the entire chain. A \code{NULL}
#'   value will initiate adaptation of \code{eps} using the dual averaging
#'   algorithm during the first \code{warmup} steps.
#' @param warmup An optional argument for how many iterations to adapt
#'   \code{eps} in the dual averaging algorithm. A value of \code{NULL}
#'   results in a default of \code{warmup=iter/2}.
#' @param adapt_delta The target acceptance rate for the dual averaging
#'   algorithm. Defaults to 80\%. NUTS does not include an accept/reject
#'   Metropolis step, so this rate can be understood as the
#'   "average acceptance probability that HMC would give to the position-momentum states explored during the final doubling iteration."
#' @param fn A function that returns the log of the posterior density.
#' @param gr A function that returns a vector of gradients of the log of
#'   the posterior density (same as \code{fn}).
#' @param covar An optional covariance matrix which can be used to improve
#'   the efficiency of sampling. The lower Cholesky decomposition of this
#'   matrix is used to transform the parameter space. If the posterior is
#'   approximately multivariate normal and \code{covar} approximates the
#'   covariance, then the transformed parameter space will be close to
#'   multivariate standard normal. In this case the algorithm will be more
#'   efficient, but there will be overhead in the matrix calculations which
#'   need to be done at each step. The default of NULL specifies to not do
#'   this transformation.
#' @param init A vector of initial parameter values.
#' @param max_treedepth Integer representing the maximum times the path
#'   length should double within an MCMC iteration. Default of 4, so 16
#'   steps. If a U-turn has not occured before this many steps the
#'   algorithm will stop and return a sample from the given tree.
#' @references \itemize{ \item{Neal, R. M. (2011). MCMC using Hamiltonian
#'   dynamics. Handbook of Markov Chain Monte Carlo.}  \item{Hoffman and
#'   Gelman (2014). The No-U-Turn sampler: Adaptively setting path lengths
#'   in Hamiltonian Monte Carlo. J. Mach. Learn. Res.  15:1593-1623.}  }
#' @return A list containing samples ('par') and algorithm details such as
#'   step size adaptation and acceptance probabilities per iteration
#'   ('sampler_params').
#' @seealso \code{\link{run_mcmc}}, \code{\link{run_mcmc.hmc}},
#'   \code{\link{run_mcmc.rwm}}
run_mcmc.nuts <- function(iter, fn, gr, init, max_treedepth=10,
                          eps=NULL, warmup=floor(iter/2),
                          adapt_delta=0.8, covar=NULL, chain=1, thin=1){
  ## If using covariance matrix and Cholesky decomposition, redefine
  ## these functions to include this transformation. The algorithm will
  ## work in the transformed space
  if(!is.null(covar)){
    fn2 <- function(theta) fn(chd %*% theta)
    gr2 <- function(theta) as.vector( t( gr(chd %*% theta) ) %*% chd )
    chd <- t(chol(covar))               # lower triangular Cholesky decomp.
    chd.inv <- solve(chd)               # inverse
    theta.cur <- chd.inv %*% init
  } else {
    fn2 <- fn; gr2 <- gr
    theta.cur <- init
  }
  sampler_params <- matrix(numeric(0), nrow=iter, ncol=6,
      dimnames=list(NULL, c("accept_stat__", "stepsize__", "treedepth__",
                            "n_leapfrog__", "divergent__", "energy__")))
  theta.out <- matrix(NA, nrow=iter, ncol=length(theta.cur))
  ## how many steps were taken at each iteration, useful for tuning
  j.results <- lp <- rep(NA, len=iter)
  useDA <- is.null(eps)               # whether to use DA algorithm
  if(useDA){
    epsvec <- Hbar <- epsbar <- rep(NA, length=warmup+1)
    eps <- epsvec[1] <- epsbar[1] <-
      .find.epsilon(theta=theta.cur, fn=fn2, gr=gr2, eps=.1, verbose=FALSE)
    mu <- log(10*eps)
    Hbar[1] <- 0; gamma <- 0.05; t0 <- 10; kappa <- 0.75
  } else {
    ## dummy values to return
    epsvec <- epsbar <- Hbar <- NULL
  }
  ## Start of MCMC chain
  time.start <- Sys.time()
  message('')
  message(paste('Starting NUTS at', time.start))
  for(m in 1:iter){
    ## Initialize this iteration from previous in case divergence at first
    ## treebuilding. If successful trajectory they are overwritten
    theta.out[m,] <- theta.minus <- theta.plus <- theta0 <- theta.cur
    lp[m] <- if(m==1) fn2(theta.cur) else lp[m-1]
    r.cur <- r.plus <- r.minus <- r0 <- rnorm(length(theta.cur),0,1)
    ## Draw a slice variable u
    u <- .sample.u(theta=theta.cur, r=r.cur, fn=fn2)
    j <- 0; n <- 1; s <- 1; divergent <- 0
    ## Track steps and divergences; updated inside .buildtree
    info <- as.environment(list(n.calls=0, divergent=0))
    while(s==1) {
      v <- sample(x=c(1,-1), size=1)
      if(v==1){
        ## move in right direction
        res <- .buildtree(theta=theta.plus, r=r.plus, u=u, v=v,
                          j=j, eps=eps, theta0=theta0, r0=r0,
                          fn=fn2, gr=gr2, info=info)
        theta.plus <- res$theta.plus
        r.plus <- res$r.plus
      } else {
        ## move in left direction
        res <- .buildtree(theta=theta.minus, r=r.minus, u=u, v=v,
                          j=j, eps=eps, theta0=theta0, r0=r0,
                          fn=fn2, gr=gr2, info=info)
        theta.minus <- res$theta.minus
        r.minus <- res$r.minus
      }
      ## test whether to accept this state
      if(res$s==1) {
        if(runif(n=1, min=0,max=1) <= res$n/n){
          theta.cur <- theta.out[m,] <- res$theta.prime
          lp[m] <- fn2(theta.cur)
          ## if(any(is.na(theta.cur))) browser()
        }
      }
      n <- n+res$n
      s <- res$s*.test.nuts(theta.plus, theta.minus, r.plus, r.minus)
      j <- j+1
      ## Stop doubling if too many or it's diverged enough
      if(j>max_treedepth) {
       ## warning("j larger than max_treedepth, skipping to next m")
        break
      }
    }
    j.results[m] <- j-1
    if(useDA){
      ## Do the adapting of eps.
      if(m <= warmup){
        alpha2 <- res$alpha/res$nalpha
        if(alpha2>1 | !is.finite(alpha2)) browser()
        Hbar[m+1] <- (1-1/(m+t0))*Hbar[m] +
          (adapt_delta-alpha2)/(m+t0)
        ## If logalpha not defined, skip this updating step and use
        ## the last one.
        if(is.nan(Hbar[m+1])) Hbar[m+1] <- abs(Hbar[m])
        logeps <- mu-sqrt(m)*Hbar[m+1]/gamma
        epsvec[m+1] <- exp(logeps)
        logepsbar <- m^(-kappa)*logeps + (1-m^(-kappa))*log(epsbar[m])
        epsbar[m+1] <- exp(logepsbar)
        eps <- epsvec[m+1]
      } else {
        eps <- epsbar[warmup]
      }
    }
    ## Save adaptation info.
    sampler_params[m,] <-
      c(res$alpha/res$nalpha, eps, j, info$n.calls, info$divergent, fn2(theta.cur))
    if(m==warmup) time.warmup <- difftime(Sys.time(), time.start, units='secs')
    .print.mcmc.progress(m, iter, warmup, chain)
  } ## end of MCMC loop
  ## Back transform parameters if covar is used
  if(!is.null(covar)) {
    theta.out <- t(apply(theta.out, 1, function(x) chd %*% x))
  }
  theta.out <- cbind(theta.out, lp)
  theta.out <- theta.out[seq(1, nrow(theta.out), by=thin),]
  sampler_params <- sampler_params[seq(1, nrow(sampler_params), by=thin),]
  ndiv <- sum(sampler_params[-(1:warmup),5])
  if(ndiv>0)
    message(paste0("There were ", ndiv, " divergent transitions after warmup"))
  message(paste0("Final acceptance ratio=", sprintf("%.2f", mean(sampler_params[-(1:warmup),1])),
                 " and target is ", adapt_delta))
  if(useDA) message(paste0("Final step size=", round(epsbar[warmup], 3),
                           "; after ", warmup, " warmup iterations"))
  time.total <- difftime(Sys.time(), time.start, units='secs')
  .print.mcmc.timing(time.warmup=time.warmup, time.total=time.total)
  return(list(par=theta.out, sampler_params=sampler_params,
              time.total=time.total, time.warmup=time.warmup,
              warmup=warmup/thin, max_treedepth=max_treedepth))
}

#' Draw a slice sample for given position and momentum variables
.sample.u <- function(theta, r, fn)
  runif(n=1, min=0, max=exp(.calculate.H(theta=theta,r=r, fn=fn)))
#' Calculate the log joint density (Hamiltonian) value for given position and
#' momentum variables.
#' @details This function currently assumes iid standard normal momentum
#' variables.
.calculate.H <- function(theta, r, fn) fn(theta)-(1/2)*sum(r^2)
#' Test whether a "U-turn" has occured in a branch of the binary tree
#' created by \ref\code{.buildtree} function. Returns TRUE if no U-turn,
#' FALSE if one occurred
.test.nuts <- function(theta.plus, theta.minus, r.plus, r.minus){
  theta.temp <- theta.plus-theta.minus
  res <- (crossprod(theta.temp,r.minus) >= 0) *
    (crossprod(theta.temp, r.plus) >= 0)
  return(res)
}

#' A recursive function that builds a leapfrog trajectory using a balanced
#' binary tree.
                                        #
#' @references This is from the No-U-Turn sampler with dual averaging
#' (algorithm 6) of Hoffman and Gelman (2014).
                                        #
#' @details The function repeatedly doubles (in a random direction) until
#' either a U-turn occurs or the trajectory becomes unstable. This is the
#' 'efficient' version that samples uniformly from the path without storing
#' it. Thus the function returns a single proposed value and not the whole
#' trajectory.
#'
.buildtree <- function(theta, r, u, v, j, eps, theta0, r0, fn, gr,
                       delta.max=1000, info = environment() ){
  if(j==0){
    ## ## Useful code for debugging. Returns entire path to global env.
    ## if(!exists('theta.trajectory'))
    ##   theta.trajectory <<- data.frame(step=0, t(theta))
    ## base case, take one step in direction v
    eps2 <- v*eps
    r <- r+(eps2/2)*gr(theta)
    theta <- theta+eps2*r
    r <- r+(eps2/2)*gr(theta)
    ## verify valid trajectory. Divergences occur if H is NaN, or drifts
    ## too from from true H.
    H <- .calculate.H(theta=theta, r=r, fn=fn)
    n <- log(u) <= H
    s <- log(u) < delta.max + H
    if(!is.finite(H) | s == 0){
     info$divergent <- 1; s <- 0
    }
    ## Acceptance ratio in log space: (Hnew-Hold)
    logalpha <- H-.calculate.H(theta=theta0, r=r0, fn=fn)
    alpha <- min(exp(logalpha),1)
    info$n.calls <- info$n.calls + 1
    ## theta.trajectory <<-
    ##   rbind(theta.trajectory, data.frame(step=tail(theta.trajectory$step,1),t(theta)))
    return(list(theta.minus=theta, theta.plus=theta, theta.prime=theta, r.minus=r,
                r.plus=r, s=s, n=n, alpha=alpha, nalpha=1))
  } else {
    ## recursion - build left and right subtrees
    xx <- .buildtree(theta=theta, r=r, u=u, v=v, j=j-1, eps=eps,
                     theta0=theta0, r0=r0, fn=fn, gr=gr, info=info)
    theta.minus <- xx$theta.minus
    theta.plus <- xx$theta.plus
    theta.prime <- xx$theta.prime
    r.minus <- xx$r.minus
    r.plus <- xx$r.plus
    alpha <- xx$alpha
    nalpha <- xx$nalpha
    s <- xx$s
    nprime <- xx$n
    ## If it didn't fail, update the above quantities
    if(s==1){
      if(v== -1){
        yy <- .buildtree(theta=theta.minus, r=r.minus, u=u, v=v,
                         j=j-1, eps=eps, theta0=theta0, r0=r0,
                         fn=fn, gr=gr, info=info)
        theta.minus <- yy$theta.minus
        r.minus <- yy$r.minus
      } else {
        yy <- .buildtree(theta=theta.plus, r=r.plus, u=u, v=v,
                         j=j-1, eps=eps, theta0=theta0, r0=r0,
                         fn=fn, gr=gr, info=info)
        theta.plus <- yy$theta.plus
        r.plus <- yy$r.plus
      }
      ## This isn't in the paper but if both slice variables failed,
      ## then you get 0/0. So I skip this test. Likewise if model
      ## throwing errors, don't keep that theta.
      nprime <- yy$n+ xx$n
      alpha <- xx$alpha+yy$alpha
      nalpha <- xx$nalpha+yy$nalpha
      if(!is.finite(nprime)) {browser();nprime <- 0}
      if(nprime!=0){
        ## choose whether to keep this theta
        if(runif(n=1, min=0, max=1) <= yy$n/nprime)
          theta.prime <- yy$theta.prime
      }
      ## check for valid proposal
      test <- .test.nuts(theta.plus=theta.plus,
                         theta.minus=theta.minus, r.plus=r.plus,
                         r.minus=r.minus)
      ## if(!test) warning(paste("U turn at j=", j))
      ## check if any of the stopping conditions were met
      s <- yy$s*test
    }
    return(list(theta.minus=theta.minus, theta.plus=theta.plus,
                theta.prime=theta.prime,
                r.minus=r.minus, r.plus=r.plus, s=s, n=nprime,
                alpha=alpha, nalpha=nalpha))
  }
}

#' Estimate a reasonable starting value for epsilon (step size) for a given
#' model, for use with Hamiltonian MCMC algorithms.
                                        #
#' This is Algorithm 4 from Hoffman and Gelman (2010) and is used in the
#' dual-averaging algorithms for both HMC and NUTS to find a reasonable
#' starting value.
#' @title Estimate step size for Hamiltonian MCMC algorithms
#' @param theta An initial parameter vector.
#' @param fn A function returning the log-likelihood (not the negative of
#' it) for a given parameter vector.
#' @param gr A function returning the gradient of the log-likelihood of a
#' model.
#' @param eps A value for espilon to initiate the algorithm. Defaults to
#' 1. If this is far too big the algorithm won't work well and an
#' alternative value can be used.
#' @return Returns the "reasonable" espilon invisible, while printing how
#' many steps to reach it.
#' @details The algorithm uses a while loop and will break after 50
#' iterations.
#'
.find.epsilon <- function(theta,  fn, gr, eps=1, verbose=TRUE){
  r <- rnorm(n=length(theta), mean=0, sd=1)
  ## Do one leapfrog step
  r.new <- r+(eps/2)*gr(theta)
  theta.new <- theta+eps*r.new
  r.new <- r.new+(eps/2)*gr(theta.new)
  H1 <- .calculate.H(theta=theta, r=r, fn=fn)
  H2 <- .calculate.H(theta=theta.new, r=r.new, fn=fn)
  a <- 2*(exp(H2)/exp(H1)>.5)-1
  ## If jumped into bad region, a can be NaN so setup algorithm to keep
  ## halving eps instead of throwing error
  if(!is.finite(a)) a <- -1
  k <- 1
  ## Similarly, keep going if there are infinite values
  while (!is.finite(H1) | !is.finite(H2) | a*H2-a*H1 > -a*log(2)) {
    eps <- (2^a)*eps
    ## Do one leapfrog step
    r.new <- r+(eps/2)*gr(theta)
    theta.new <- theta+eps*r.new
    r.new <- r.new+(eps/2)*gr(theta.new)
    H2 <- .calculate.H(theta=theta.new, r=r.new, fn=fn)
    k <- k+1
    if(k>50) {
      stop("More than 50 iterations to find reasonable eps. Model is likely misspecified or some other issue.")
    }
  }
  if(verbose) message(paste("Reasonable epsilon=", eps, "found after", k, "steps"))
  return(invisible(eps))
}

#' Print MCMC progress to console.
#'
#' @param iteration The iteration of the MCMC chain.
#' @param iter The total iterations.
#' @param warmup The number of warmup iterations.
#' @param chain The chain being run (bookkeeping only).
#' @return Nothing. Prints to message to console.
#'
#' @details This function was modeled after the functionality provided by
# the R package \link{rstan}.
.print.mcmc.progress <- function(iteration, iter, warmup, chain){
  i <- iteration
  refresh <- max(10, floor(iter/10))
  if(i==1 | i==iter | i %% refresh ==0){
    i.width <- formatC(i, width=nchar(iter))
    out <- paste0('Chain ',chain,', Iteration: ', i.width , "/", iter, " [",
                  formatC(floor(100*(i/iter)), width=3), "%]",
                  ifelse(i <= warmup, " (Warmup)", " (Sampling)"))
    message(out)
  }
}

#' Print MCMC timing to console
#' @param time.warmup Time of warmup in seconds.
#' @param time.total Time of total in seconds.
#' @return Nothing. Prints message to console.
#'
#' @details This function was modeled after the functionality provided by
#'   the R package \link{rstan}.
.print.mcmc.timing <- function(time.warmup, time.total){
  x <- ' Elapsed Time: '
  message(paste0(x, sprintf("%.1f", time.warmup), ' seconds (Warmup)'))
  message(paste0(x, sprintf("%.1f", time.total-time.warmup), ' seconds (Sampling)'))
  message(paste0(x, sprintf("%.1f", time.total), ' seconds (Total)'))
}

#' Convert TMB output from \link{\code{run_mcmc}} into a \code{shinystan}
#' object.
#'
#' @details The shinystan packages provides several conversion functions
#'   for objects of different types, such as stanfit classes (Stan ouput)
#'   and simple arrays. For the latter, option NUTS information, such as
#'   \code{sampler_params} can be passed. This function essentially extends
#'   the functionality of \code{as.shinystan} to work specifically with TMB
#'   MCMC lists. The user can thus explore their TMB model with
#'   \code{launch_shinystan(as.shinystan.tmb(tmb.fit))} in the same way
#'   that Stan models are examined.
#' @param tmb.fit Output list from \link{\code{run_mcmc}} for any of the
#' three algorithms.
#' @seealso launch_shinystan_tmb
#' @return An S4 object of class shinystan. Depending on the algorithm
#'   used, this list will have slight differences.
as.shinystan.tmb <- function(tmb.fit){
  if(tmb.fit$algorithm=="NUTS"){
    sso <- with(tmb.fit, as.shinystan(samples, burnin=warmup, max_treedepth=max_treedepth,
             sampler_params=sampler_params, algorithm='NUTS', model_name=model))
  } else if(tmb.fit$algorithm=="HMC"){
    sso <- with(tmb.fit, as.shinystan(samples, burnin=warmup,
             sampler_params=sampler_params, algorithm='HMC', model_name=model))
  } else {
    sso <- with(tmb.fit, as.shinystan(samples, burnin=warmup,
             algorithm='RWM', model_name=model))
  }
  return(sso)
}

#' A high level wrapper to launch shinystan for a TMB MCMC list object.
#'
#' @details This function simply calls
#'   \code{launch_shinystan(as.shinystan.tmb(tmb.fit))}.
launch_shinystan_tmb <- function(tmb.fit){
  launch_shinystan(as.shinystan.tmb(tmb.fit))
}

#' Extract posterior samples from a TMB MCMC fit list.
#'
#' @param fit.tmb A list returned by \code{\link{run_mcmc}}.
#' @param inc_warmup Whether to extract the warmup samples or not
#'   (default). Warmup samples should never be used for inference, but may
#'   be useful for diagnostics.
#' @return An invisible data.frame containing samples (rows) of each
#'   parameter (columns). If multiple chains exist they will be rbinded
#'   together.
extract_samples <- function(fit.tmb, inc_warmup=FALSE){
  x <- fit.tmb$samples
  if(!is.array(x)) stop("fit.tmb$samples is not an array -- valid TMB output?")
  ind <- if(inc_warmup) 1:dim(x)[1] else -(1:fit.tmb$warmup)
  y <- do.call(rbind, lapply(1:dim(x)[2], function(i) x[ind, i, -dim(x)[3]]))
  return(invisible(as.data.frame(y)))
}
