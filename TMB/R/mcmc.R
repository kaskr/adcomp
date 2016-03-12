## Copyright (C) 2015 Cole Monnahan
## License: GPL-2

#' [BETA VERSION] Draw samples from the posterior of a TMB model using a
#' specified MCMC algorithm.
#'
#' @details This function is a top-level wrapper designed specifically to
#' work with TMB models. There are several MCMC algorithms available for
#' use. The user is responsible for specifying the model properly (priors,
#' starting values, desired parameters fixed, etc.), as well as assessing
#' the convergence and validity of the resulting samples (e.g., through the
#' \code{coda} package) before making inference.
#' @title MCMC sampling of TMB models
#' @author Cole Monnahan
#' @param obj A TMB model object.
#' @param nsim The number of (dependent) samples to draw.
#' @param params.init The initial parameter vector. The default of NULL
#' signifies to use the starting values present in the model
#' (i.e., \code{obj$par}).
#' @param covar An optional covariance matrix which can be used to improve
#' the efficiency of sampling. The lower Cholesky decomposition of this
#' matrix is used to transform the parameter space. If the posterior is
#' approximately multivariate normal and \code{covar} approximates the
#' covariance, then the transformed parameter space will be close to
#' multivariate standard normal. In this case the algorithm will be more
#' efficient, but there will be overhead in the matrix calculations which
#' need to be done at each step. The default of NULL specifies to not do
#' this transformation.
#' @param algorithm A string specifiying an algorithm. Currently supported
#' are: \itemize{
#' \item{"RWM"}{the random walk Metropolis sampler}
#' \item{"HMC"}{the Hamiltonian sampler (see Neal 2011)}
#' \item{"NUTS"}{the No-U-Turn sampler (see Hoffman and Gelman 2014)}
#' }
#' These algorithms require different arguments; see their help files for more information.
#' @param diagnostic Whether to return diagnostic information about
#' chain. See individual algorithm for more information.
#' @param ... Further arguments to be passed to the algorithm. See help
#' files for the samplers for further arguments.
#' @return If \code{diagnostic} is \code{FALSE}, returns a data frame with
#' posterior samples. Otherwise it returns a list containing the samples
#' and properties of the sampler useful for diagnosing behavior and
#' efficiency.
#' @example inst/examples/mcmc_examples.R
#' @seealso \code{\link{run_mcmc.hmc}}, \code{\link{run_mcmc.nuts}}, \code{\link{run_mcmc.rwm}}
run_mcmc <- function(obj, nsim, algorithm, params.init=NULL, covar=NULL, diagnostic=FALSE, ...){
    ## Initialization for all algorithms
    algorithm <- match.arg(algorithm, choices=c("HMC", "NUTS", "RWM"))
    fn <- function(x) {
        z <- -obj$fn(x)
        if(is.nan(z)){
            warning(paste("replacing NaN w/ Inf at:", paste(x, collapse=" ")))
            z <- Inf
        }
        return(z)
    }
    gr <- function(x) {
        z <- -as.vector(obj$gr(x))
        if(any(is.nan(z))){
            warning(paste("NaN gradient at:", paste(x, collapse=" ")))
            z <- rep(0, length(x))
           }
        return(z)
    }
    obj$env$beSilent()                  # silence console output
    ## argument checking
    if(is.null(params.init)){
        params.init <- obj$par
    } else if(length(params.init) != length(obj$par)){
        stop("params.init is wrong length")
    }
    ## Select and run the chain.
    if(algorithm=="HMC")
        time <- system.time(mcmc.out <-
            run_mcmc.hmc(nsim=nsim, fn=fn, gr=gr, params.init=params.init, covar=covar,
                     diagnostic=diagnostic, ...))
    else if(algorithm=="NUTS")
        time <- system.time(mcmc.out <-
            run_mcmc.nuts(nsim=nsim, fn=fn, gr=gr, params.init=params.init, covar=covar,
                      diagnostic=diagnostic, ...))
    else if(algorithm=="RWM")
        time <- system.time(mcmc.out <-
            run_mcmc.rwm(nsim=nsim, fn=fn, params.init=params.init, covar=covar,
                      diagnostic=diagnostic, ...))
    ## Clean up returned output, a matrix if diag is FALSE, otherwise a list
    if(!diagnostic){
        mcmc.out <- as.data.frame(mcmc.out)
        names(mcmc.out) <- names(obj$par)
    } else {
        mcmc.out$time <- as.numeric(time[3])        # grab the elapsed time
        mcmc.out$par <- as.data.frame(mcmc.out$par)
        names(mcmc.out$par) <- names(obj$par)
    }
    return(invisible(mcmc.out))
}


#' [BETA VERSION] Draw MCMC samples from a model posterior using a
#' Random Walk Metropolis (RWM) sampler.
#'
#' @param nsim The number of samples to return.
#' @param fn A function that returns the log of the posterior density.
#' @param params.init A vector of initial parameter values.
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
#' \code{nsim} samples from the posterior. Otherwise returns a list
#' containing samples ('par'), proposed samples ('par.proposed'), vector of
#' which proposals were accepted ('accepted'), and the total function calls
#' ('n.calls'), which for this algorithm is \code{nsim}
#' @seealso \code{\link{run_mcmc}}, \code{\link{run_mcmc.nuts}}, \code{\link{run_mcmc.hmc}}
run_mcmc.rwm <- function(nsim, fn, params.init, alpha=1, covar=NULL, diagnostic=FALSE){
    accepted <- rep(0, length=nsim)
    n.params <- length(params.init)
    theta.out <- matrix(NA, nrow=nsim, ncol=n.params)
    if(diagnostic) theta.proposed <- theta.out
    ## If using covariance matrix and Cholesky decomposition, redefine
    ## these functions to include this transformation. The algorithm will
    ## work in the transformed space.
    if(!is.null(covar)){
        fn2 <- function(theta) fn(chd %*% theta)
        chd <- t(chol(covar))               # lower triangular Cholesky decomp.
        chd.inv <- solve(chd)               # inverse
        theta.cur <- chd.inv %*% params.init
    } else {
        fn2 <- fn
        theta.cur <- params.init
    }
    fn.cur <- fn2(theta.cur)
    for(m in 1:nsim){
        ## generate proposal
        theta.new <- theta.cur + alpha*rnorm(n=n.params, mean=0, sd=1)
        fn.new <- fn2(theta.new)
        if(diagnostic) theta.proposed[m,] <- theta.new
        if(log(runif(1))< fn.new-fn.cur){
            ## accept
            accepted[m] <- 1
            theta.cur <- theta.out[m,] <- theta.new
            fn.cur <- fn.new
        } else {
            ## do not accept
            theta.out[m,] <- theta.cur
        }
    } # end of MCMC loop
    ## Back transform parameters if covar is used
    if(!is.null(covar)) {
        theta.out <- t(apply(theta.out, 1, function(x) chd %*% x))
        if(diagnostic)
            theta.proposed <- t(apply(theta.proposed, 1, function(x) chd %*% x))
    }
    message(paste("Acceptance rate = ", round(mean(accepted),1)))
    if(diagnostic){
        theta.out <- list(par=theta.out, accepted=accepted,
                         acceptance=mean(accepted), n.calls=nsim,
                         par.proposed=theta.proposed)
    }
    return(theta.out)
}

#' [BETA VERSION] Draw MCMC samples from a model posterior using a
#' Hamiltonian sampler.
#' @details This function implements algorithm 5 of Hoffman and Gelman
#' (2014), which includes adaptive step sizes (\code{eps}) via an algorithm
#' called dual averaging. In theory neither the step length nor step size
#' needs to be input by the user to obtain efficient sampling from the
#' posterior.
#' @param nsim The number of samples to return.
#' @param L The number of leapfrog steps to take. The NUTS algorithm does
#' not require this as an input. If \code{L=1} this function will perform
#' Langevin sampling. In some contexts \code{L} can roughly be thought of
#' as a thinning rate.
#' @param eps The length of the leapfrog steps. If a numeric value is
#' passed, it will be used throughout the entire chain. A \code{NULL}
#' value will initiate adaptation of \code{eps} using the dual averaging
#' algorithm during the first \code{Madapt} steps.
#' @param Madapt An optional argument for how many iterations to adapt
#' \code{eps} in the dual averaging algorithm. A value of \code{NULL}
#' results in a default of \code{Madapt=nsim/2}.
#' @param delta The target acceptance rate if using apative
#' \code{eps}. Defaults to 50\%.
#' @param fn A function that returns the log of the posterior density.
#' @param gr A function that returns a vector of gradients of the log of
#' the posterior density (same as \code{fn}).
#' @param params.init A vector of initial parameter values.
#' @param covar An optional covariance matrix which can be used to improve
#' the efficiency of sampling. The lower Cholesky decomposition of this
#' matrix is used to transform the parameter space. If the posterior is
#' approximately multivariate normal and \code{covar} approximates the
#' covariance, then the transformed parameter space will be close to
#' multivariate standard normal. In this case the algorithm will be more
#' efficient, but there will be overhead in the matrix calculations which
#' need to be done at each step. The default of NULL specifies to not do
#' this transformation.
#' @param diagnostic Whether to return a list of diagnostic metrics about
#' the chain. Useful for assessing efficiency and tuning chain.
#' @references
#' \itemize{
#' \item{Neal, R. M. (2011). MCMC using Hamiltonian dynamics. Handbook
#' of Markov Chain Monte Carlo.}
#' \item{Hoffman and Gelman (2014). The No-U-Turn sampler: Adaptively
#' setting path lengths in Hamiltonian Monte Carlo. J. Mach. Learn. Res.
#' 15:1593-1623.}
#' }
#' @seealso \code{\link{run_mcmc}}, \code{\link{run_mcmc.nuts}}, \code{\link{run_mcmc.rwm}}
#' @return If \code{diagnostic} is FALSE (default), returns a matrix of
#' \code{nsim} samples from the posterior. Otherwise returns a list
#' containing samples ('par'), proposed samples ('par.proposed'), vector of
#' which were accepted ('accepted'), and the total function and gradient
#' calls ('n.calls'), which for this algorithm is \code{nsim*(L+2)}.
run_mcmc.hmc <- function(nsim, fn, gr, params.init, L, eps=NULL, covar=NULL,
                     delta=0.5, Madapt=NULL, diagnostic=FALSE){
    ## If using covariance matrix and Cholesky decomposition, redefine
    ## these functions to include this transformation. The algorithm will
    ## work in the transformed space
    if(!is.null(covar)){
        fn2 <- function(theta) fn(chd %*% theta)
        gr2 <- function(theta) as.vector( t( gr(chd %*% theta) ) %*% chd )
        chd <- t(chol(covar))               # lower triangular Cholesky decomp.
        chd.inv <- solve(chd)               # inverse
        theta.cur <- chd.inv %*% params.init
    } else {
        fn2 <- fn; gr2 <- gr
        theta.cur <- params.init
    }
    accepted <- rep(NA, nsim)
    theta.out <- matrix(NA, nrow=nsim, ncol=length(params.init))
    if(diagnostic) theta.proposed <- theta.out
    ## A NULL value for eps signifies to use the dual averaging algorithm
    useDA <- is.null(eps)
    if(useDA){
        if(is.null(Madapt)){
            message("MCMC HMC: Madapt not specified, defaulting to half of nsim")
            Madapt <- floor(nsim/2)
        }
        ## Initialize the dual-averaging algorithm.
        message(paste("MCMC HMC: No eps given so using dual averaging during first", Madapt, "steps."))
        epsvec <- Hbar <- epsbar <- rep(NA, length=Madapt+1)
        eps <- epsvec[1] <- epsbar[1] <-
            .find.epsilon(theta=theta.cur, fn=fn2, gr=gr2, eps=.1, verbose=FALSE)
        mu <- log(10*eps)
        Hbar[1] <- 0; gamma <- 0.05; t0 <- 10; kappa <- 0.75
    } else {
        ## dummy values to return
        epsvec <- epsbar <- Hbar <- NULL
    }
    ## Start of MCMC chain
    eps0 <- eps                         # eps0 is average eps
    for(m in 1:nsim){
        ## Randomize eps if not doing adapation below
        if(!useDA) eps <- eps0*runif(1,.9,1.1)
        r.cur <- r.new <- rnorm(length(params.init),0,1)
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
        }
        ## half step for momentum at the end
        r.new <- r.new+eps*gr2(theta.new)/2
        ## negate r to make proposal symmetric
        r.new <- -r.new
        if(diagnostic) theta.proposed[m,] <- theta.new
        logalpha <- -fn2(theta.cur)+fn2(theta.new)+ sum(r.cur^2)/2-sum(r.new^2)/2
        if(is.finite(logalpha) & log(runif(1)) < logalpha){
            ## accept the proposed state
            theta.cur <- theta.new
            accepted[m] <- TRUE
        } else {
            ## otherwise reject it and stay there
            accepted[m] <- FALSE
        }
        theta.out[m,] <- theta.cur
        if(useDA){
            ## Do the adapting of eps.
            if(m <= Madapt){
                Hbar[m+1] <-
                    (1-1/(m+t0))*Hbar[m] +
                        (delta-min(1,exp(logalpha)))/(m+t0)
                ## If logalpha not defined, skip this updating step and use
                ## the last one.
                if(is.nan(Hbar[m+1])) Hbar[m+1] <- abs(Hbar[m])
                logeps <- mu-sqrt(m)*Hbar[m+1]/gamma
                epsvec[m+1] <- exp(logeps)
                logepsbar <- m^(-kappa)*logeps + (1-m^(-kappa))*log(epsbar[m])
                epsbar[m+1] <- exp(logepsbar)
                eps <- epsvec[m+1]
            } else {
                eps <- epsbar[Madapt]*runif(1,.9,1.1)
            }
        }
    } ## end of MCMC loop
    ## Back transform parameters if covar is used
    if(!is.null(covar)) {
        theta.out <- t(apply(theta.out, 1, function(x) chd %*% x))
        if(diagnostic)
            theta.proposed <- t(apply(theta.proposed, 1, function(x) chd %*% x))
    }
    message(paste("MCMC HMC: Acceptance rate = ", round(mean(accepted),1)))
    if(useDA) message(paste("MCMC HMC: Dual averaging final average eps =", round(epsbar[Madapt], 3)))
    if(diagnostic){
        return(list(par=theta.out, par.proposed=theta.proposed, accepted=accepted,
                    n.calls=nsim*(L+2), epsvec=epsvec, epsbar=epsbar, Hbar=Hbar))
    } else {
        return(theta.out)
    }
}


#' [BETA VERSION] Draw MCMC samples from a model posterior using the
#' No-U-Turn (NUTS) sampler with dual averaging.
#'
#' @details This function implements algorithm 6 of Hoffman and Gelman
#' (2014), which includes adaptive step sizes (\code{eps}) via an algorithm
#' called dual averaging. In theory neither the step length nor step size
#' needs to be input by the user to obtain efficient sampling from the
#' posterior.
#' @param nsim The number of samples to return.
#' @param eps The length of the leapfrog steps. If a numeric value is
#' passed, it will be used throughout the entire chain. A \code{NULL}
#' value will initiate adaptation of \code{eps} using the dual averaging
#' algorithm during the first \code{Madapt} steps.
#' @param Madapt An optional argument for how many iterations to adapt
#' \code{eps} in the dual averaging algorithm. A value of \code{NULL}
#' results in a default of \code{Madapt=nsim/2}.
#' @param delta The target acceptance rate for the dual averaging
#' algorithm. Defaults to 50\%. NUTS does not include an accept/reject
#' Metropolis step, so this rate can be understood as the "average acceptance
#' probability that HMC would give to the position-momentum states explored
#' during the final doubling iteration."
#' @param fn A function that returns the log of the posterior density.
#' @param gr A function that returns a vector of gradients of the log of
#' the posterior density (same as \code{fn}).
#' @param covar An optional covariance matrix which can be used to improve
#' the efficiency of sampling. The lower Cholesky decomposition of this
#' matrix is used to transform the parameter space. If the posterior is
#' approximately multivariate normal and \code{covar} approximates the
#' covariance, then the transformed parameter space will be close to
#' multivariate standard normal. In this case the algorithm will be more
#' efficient, but there will be overhead in the matrix calculations which
#' need to be done at each step. The default of NULL specifies to not do
#' this transformation.
#' @param params.init A vector of initial parameter values.
#' @param diagnostic Whether to return a list of diagnostic metrics about
#' the chain. Useful for assessing efficiency and tuning chain.
#' @param max_doublings Integer representing the maximum times the path
#' length should double within an MCMC iteration. Default of 4, so 16
#' steps. If a U-turn has not occured before this many steps the algorithm
#' will stop and return a sample from the given tree.
#' @references
#'  \itemize{
#' \item{Neal, R. M. (2011). MCMC using Hamiltonian dynamics. Handbook
#' of Markov Chain Monte Carlo.}
#' \item{Hoffman and Gelman (2014). The No-U-Turn sampler: Adaptively
#' setting path lengths in Hamiltonian Monte Carlo. J. Mach. Learn. Res.
#' 15:1593-1623.}
#' }
#' @return If \code{diagnostic} is FALSE (default), returns a matrix
#' of \code{nsim} samples from the posterior. Otherwise returns a list
#' containing samples ('par'), vector of steps taken at each iteration
#' ('steps.taken'), and the total function and gradient calls
#' ('n.calls'), which in the case of NUTS is dynamic, and finally the
#' average \code{eps} ('epsbar') from the dual averaging algorithm if
#' used (otherwise NULL).
#' @seealso \code{\link{run_mcmc}}, \code{\link{run_mcmc.hmc}}, \code{\link{run_mcmc.rwm}}
run_mcmc.nuts <- function(nsim, fn, gr, params.init, max_doublings=4, eps=NULL, Madapt=NULL,
                      delta=0.5, covar=NULL, diagnostic=FALSE){
    ## If using covariance matrix and Cholesky decomposition, redefine
    ## these functions to include this transformation. The algorithm will
    ## work in the transformed space
    if(!is.null(covar)){
        fn2 <- function(theta) fn(chd %*% theta)
        gr2 <- function(theta) as.vector( t( gr(chd %*% theta) ) %*% chd )
        chd <- t(chol(covar))               # lower triangular Cholesky decomp.
        chd.inv <- solve(chd)               # inverse
        theta.cur <- chd.inv %*% params.init
    } else {
        fn2 <- fn; gr2 <- gr
        theta.cur <- params.init
    }
    theta.out <- matrix(NA, nrow=nsim, ncol=length(theta.cur))
    ## how many steps were taken at each iteration, useful for tuning
    j.results <- rep(NA, len=nsim)
    ## count the model calls; updated inside .buildtree. Some subtrees
    ## wont finish due to exit conditions so this is dynamic and not a
    ## simple formula like with HMC.
    info <- as.environment( list(n.calls = 0) )
    useDA <- is.null(eps)               # whether to use DA algorithm
    if(useDA){
        if(is.null(Madapt)){
            message("MCMC NUTS: Madapt not specified, defaulting to half of nsim")
            Madapt <- floor(nsim/2)
        }
        ## Initialize the dual-averaging algorithm.
        message(paste("MCMC NUTS: No eps given so using dual averaging during first", Madapt, "steps."))
        epsvec <- Hbar <- epsbar <- rep(NA, length=Madapt+1)
        eps <- epsvec[1] <- epsbar[1] <-
            .find.epsilon(theta=theta.cur, fn=fn2, gr=gr2, eps=.1, verbose=FALSE)
        mu <- log(10*eps)
        Hbar[1] <- 0; gamma <- 0.05; t0 <- 10; kappa <- 0.75
    } else {
        ## dummy values to return
        epsvec <- epsbar <- Hbar <- NULL
    }
    ## Start of MCMC chain
    for(m in 1:nsim){
        ## initialize
        theta.out[m,] <- theta.minus <- theta.plus <- theta0 <- theta.cur
        r.cur <- r.plus <- r.minus <- r0 <- rnorm(length(theta.cur),0,1)
        ## Draw a slice variable u
        u <- .sample.u(theta=theta.cur, r=r.cur, fn=fn2)
        j <- 0; n <- 1; s <- 1
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
            if(is.na(res$s) | is.nan(res$s))  res$s <- 0
            if(res$s==1) {
                if(runif(n=1, min=0,max=1) <= res$n/n){
                    theta.cur <- res$theta.prime
                    theta.out[m,] <- res$theta.prime
                }
            }
            n <- n+res$n
            s <- res$s*.test.nuts(theta.plus, theta.minus, r.plus, r.minus)
            ## Stop trajectory if there are any problems, probably happens
            ## when jumping way too far into the tails and the model isn't
            ## defined
            if(is.na(s) | is.nan(s))  s <- 0
            j <- j+1
            ## Stop doubling if too many or it's diverged enough
            if(j>max_doublings & s) {
                ## warning("j larger than max_doublings, skipping to next m")
                break
            }
        }
        j.results[m] <- j-1
        if(useDA){
            ## Do the adapting of eps.
            if(m <= Madapt){
                Hbar[m+1] <- (1-1/(m+t0))*Hbar[m] +
                    (delta-res$alpha/res$nalpha)/(m+t0)
                ## If logalpha not defined, skip this updating step and use
                ## the last one.
                if(is.nan(Hbar[m+1])) Hbar[m+1] <- abs(Hbar[m])

                logeps <- mu-sqrt(m)*Hbar[m+1]/gamma
                epsvec[m+1] <- exp(logeps)
                logepsbar <- m^(-kappa)*logeps + (1-m^(-kappa))*log(epsbar[m])
                epsbar[m+1] <- exp(logepsbar)
                eps <- epsvec[m+1]
            } else {
                eps <- epsbar[Madapt]*runif(1,.9,1.1)
            }
        }
    } ## end of MCMC loop
    ## Back transform parameters if covar is used
    if(!is.null(covar)) {
        theta.out <- t(apply(theta.out, 1, function(x) chd %*% x))
    }
    j.stats <- 2^(c(min(j.results), median(j.results), max(j.results)))
    if(useDA)
        message(paste("MCMC NUTS: Dual averaging final average eps =", round(epsbar[Madapt], 3)))
    message(paste0("MCMC NUTS: Approximate leapfrog steps(min, median, max)=(",
                   paste(j.stats, collapse=","), ")"))
    if(diagnostic){
        return(list(par=theta.out, steps.taken= 2^j.results,
                    n.calls=info$n.calls, epsvec=epsvec, epsbar=epsbar, Hbar=Hbar))
    } else {
        return(theta.out)
    }
}

# Draw a slice sample for given position and momentum variables
.sample.u <- function(theta, r, fn)
    runif(n=1, min=0, max=exp(.calculate.H(theta=theta,r=r, fn=fn)))
# Calculate the Hamiltonian value for position and momentum variables.
#
# @details This function currently assumes iid standard normal momentum
# variables.
.calculate.H <- function(theta, r, fn) fn(theta)-(1/2)*sum(r^2)
# Test whether a "U-turn" has occured in a branch of the binary tree
# created by \ref\code{.buildtree} function.
.test.nuts <- function(theta.plus, theta.minus, r.plus, r.minus){
    theta.temp <- (theta.plus-theta.minus)
   as.numeric( theta.temp %*% r.minus >= 0 | theta.temp %*% r.plus >= 0)
}

# A recursive function that builds a leapfrog trajectory using a balanced
# binary tree.
#
# @references This is from the No-U-Turn sampler with dual averaging
# (algorithm 6) of Hoffman and Gelman (2014).
#
# @details The function repeatedly doubles (in a random direction) until
# either a U-turn occurs or the trajectory becomes unstable. This is the
# 'efficient' version that samples uniformly from the path without storing
# it. Thus the function returns a single proposed value and not the whole
# trajectory.
.buildtree <- function(theta, r, u, v, j, eps, theta0, r0, fn, gr,
                         delta.max=1000, info = environment() ){
    if(j==0){
        ## base case, take one step in direction v
        eps <- v*eps
        r <- r+(eps/2)*gr(theta)
        theta <- theta+eps*r
        r <- r+(eps/2)*gr(theta)
        ## verify valid trajectory
        H <- .calculate.H(theta=theta, r=r, fn=fn)
        s <- H-log(u) + delta.max > 0
        if(is.na(s) | is.nan(s)) s <- 0
        n <- log(u) <= H
        ## ## Useful code for debugging. Returns entire path to global env.
        ## if(!exists('theta.trajectory'))
        ##     theta.trajectory <<- theta
        ## else
        ##     theta.trajectory <<- rbind(theta.trajectory, theta)
        temp <- .calculate.H(theta=theta, r=r, fn=fn)-
            .calculate.H(theta=theta0, r=r0, fn=fn)
        alpha <- min(exp(temp),1)
        info$n.calls <- info$n.calls + 5
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
        if(is.na(s) | is.nan(s)) s <- 0
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
            if(!is.finite(nprime)) nprime <- 0
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
            s <- xx$s*yy$s*test
        }
        return(list(theta.minus=theta.minus, theta.plus=theta.plus,
                    theta.prime=theta.prime,
                    r.minus=r.minus, r.plus=r.plus, s=s, n=nprime,
                    alpha=alpha, nalpha=1))
    }
}

# Estimate a reasonable starting value for epsilon (step size) for a given
# model, for use with Hamiltonian MCMC algorithms.
#
# This is Algorithm 4 from Hoffman and Gelman (2010) and is used in the
# dual-averaging algorithms for both HMC and NUTS to find a reasonable
# starting value.
# @title Estimate step size for Hamiltonian MCMC algorithms
# @param theta An initial parameter vector.
# @param fn A function returning the log-likelihood (not the negative of
# it) for a given parameter vector.
# @param gr A function returning the gradient of the log-likelihood of a
# model.
# @param eps A value for espilon to initiate the algorithm. Defaults to
# 1. If this is far too big the algorithm won't work well and an
# alternative value can be used.
# @return Returns the "reasonable" espilon invisible, while printing how
# many steps to reach it.
# @details The algorithm uses a while loop and will break after 50
# iterations.
#
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
