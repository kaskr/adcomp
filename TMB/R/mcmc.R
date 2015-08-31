#' [BETA VERSION] Draw samples from the posterior of a TMB model using a
#' specified MCMC algorithm.
#'
#' @details The user is responsible for specifying the model properly
#' (priors, starting values, desired parameters fixed, etc.).
#' @title MCMC sampling of TMB models
#' @author Cole Monnahan
#' @param obj A TMB model object.
#' @param nsim The number of (dependent) samples to draw.
#' @param params.init The initial parameter vector. The default of NULL
#' signifies to use the starting values present in the model
#' (i.e., \code{obj$par}).
#' @param algorithm A string specifiying an algorithm. Currently supported
#' are "HMC" for Hamiltonian sampler and "NUTS" for the No-U-Turn
#' sampler.
#' @param diagnostic Whether to return diagnostic information about
#' chain. See individual algorithm for more information.
#' @param ... Further arguments to be passed to the algorithm. See help
#' files for the samplers for further arguments.
#' @example inst/examples/mcmc_examples.R
run_mcmc <- function(obj, nsim, algorithm, params.init=NULL, diagnostic=FALSE, ...){
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
               warning(paste("NaN at:", paste(x, collapse=" ")))
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
            mcmc.hmc(nsim=nsim, fn=fn, gr=gr, params.init=params.init,
                     diagnostic=diagnostic, ...))
    else if(algorithm=="NUTS")
        time <- system.time(mcmc.out <-
            mcmc.nuts(nsim=nsim, fn=fn, gr=gr, params.init=params.init,
                      diagnostic=diagnostic, ...))
    else if(algorithm=="RWM")
        time <- system.time(mcmc.out <-
            mcmc.rwm(nsim=nsim, fn=fn, params.init=params.init,
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
#' Random Walk Metropolis sampler.
#'
#' @param nsim The number of samples to return.
#' @param fn A function that returns the log of the posterior density.
#' @param params.init A vector of initial parameter values.
#' @param diagnostic Whether to return a list of diagnostic metrics about
#' the chain. Useful for assessing efficiency and tuning chain.
#' @details This
#' @param covar A covariance matrix to be used in generating multivariate
#' normal proposals. A value of NULL (default) indicates to use iid
#' standard normal.
#' @param alpha The amount to scale the proposal, i.e,
#' Xnew=Xcur+alpha*Xproposed where Xproposed is generated from a mean-zero
#' multivariate normal. Varying \code{alpha} varies the acceptance rate.
#' @return If \code{diagnostic} is FALSE (default), returns a matrix of
#' \code{nsim} samples from the posterior. Otherwise returns a list
#' containing samples ('par'), proposed samples ('par.proposed'), vector of
#' which were accepted ('accepted'), and the total function calls
#' ('n.calls'), which for this algorithm is \code{nsim}
mcmc.rwm <- function(nsim, fn, params.init, alpha=1, covar=NULL, diagnostic=FALSE){
    accepted <- rep(0, length=nsim)
    n.params <- length(params.init)
    theta.out <- matrix(NA, nrow=nsim, ncol=n.params)
    if(diagnostic) theta.proposed <- theta.out
    theta.cur <- theta.out[1,] <- params.init
    if(is.null(covar)) f <- function() rnorm(n=n.params, 0, 1)
        else f <-  function() mvtnorm::rmvnorm(n=1, mean=rep(0, n.params), sigma=covar)
    for(m in 2:nsim){
        ## generate proposal
        theta.new <- theta.cur + alpha*f()
        if(diagnostic) theta.proposed[m,] <- theta.new
        if(log(runif(1))< fn(theta.new)-fn(theta.cur)){
            ## accept
            accepted[m] <- 1
            theta.cur <- theta.out[m,] <- theta.new
        } else {
            ## do not accept
            theta.out[m,] <- theta.cur
        }
    }
    if(diagnostic){
        theta.out <- list(par=theta.out, accepted=accepted,
                         acceptance=mean(accepted), n.calls=nsim,
                         par.proposed=theta.proposed)
    }
    return(theta.out)
}



#' [BETA VERSION] Draw MCMC samples from a model posterior using a Hamiltonian sampler.
#'
#' @param nsim The number of samples to return.
#' @param L The number of leapfrog steps to take. The NUTS algorithm does
#' not require this as an input. If L is 1 this function will perform
#' Langevin sampling.
#' @param eps The length of the leapfrog steps.
#' @param fn A function that returns the log of the posterior density. This
#' function should not return the negative log, following the notation of
#' Hoffman and Gelman (2014).
#' @param gr A function that returns a vector of gradients of the log of
#' the posterior density (same as with \code{fn}).
#' @param params.init A vector of initial parameter values.
#' @param diagnostic Whether to return a list of diagnostic metrics about
#' the chain. Useful for assessing efficiency and tuning chain.
#' @references Neal, R. M. 2011. MCMC using Hamiltonian dynamics.
#' @return If \code{diagnostic} is FALSE (default), returns a matrix of
#' \code{nsim} samples from the posterior. Otherwise returns a list
#' containing samples ('par'), proposed samples ('par.proposed'), vector of
#' which were accepted ('accepted'), and the total function and gradient
#' calls ('n.calls'), which for this algorithm is \code{nsim}*(\code{L}+2)
mcmc.hmc <- function(nsim, L, eps, fn, gr, params.init, covar=NULL,
                      diagnostic=FALSE){
    ## If using covariance matrix and Cholesky decomposition, redefine
    ## these functions to include this transformation. The algorithm will
    ## work in the transformed space
    if(!is.null(covar)){
        fn2 <- function(theta) fn(chd %*% theta)
        gr2 <- function(theta) gr(chd %*% theta)
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
    eps0 <- eps
    for(m in 1:nsim){
        eps <- eps0*runif(1,.9,1.1)
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
        alpha <- -fn2(theta.cur)+fn2(theta.new)+ sum(r.cur^2)/2-sum(r.new^2)/2
        if(is.finite(alpha) & log(runif(1)) < alpha){
            ## accept the proposed state
            theta.cur <- theta.new
            accepted[m] <- TRUE
        } else {
            ## otherwise reject it and stay there
            accepted[m] <- FALSE
        }
        theta.out[m,] <- theta.cur
    }
    ## Back transform parameters if covar is used
    if(!is.null(covar)) {
        theta.out <- t(apply(theta.out, 1, function(x) chd %*% x))
        if(diagnostic)
            theta.proposed <- t(apply(theta.proposed, 1, function(x) chd %*% x))
    }
    message(paste("Acceptance rate = ", round(mean(accepted),1)))
    if(diagnostic){
        return(list(par=theta.out, par.proposed=theta.proposed, accepted=accepted,
                    n.calls=nsim*(L+2)))
    } else {
        return(theta.out)
    }
}


#' [BETA VERSION] Draw MCMC samples from a model posterior using the
#' No-U-Turn (NUTS) sampler with dual averaging from Hoffman and Gelman
#' (2014).
#'
#' @details This is the 'efficient' NUTS algorithm but does not use dual
#' averaging to tune \code{eps}, so it must be specified.
#' @param nsim The number of samples to return.
#' @param eps The length of the leapfrog steps. If NULL is passed (the
#' default) then dual averaging is used during the first \code{Madapt}
#' steps.
#' @param delta The target acceptance rate for the dual averaging
#' algorithm. Must be specified if \code{eps} is NULL. The rate can be
#' "understood as the average acceptance probability that HMC would give to
#' the position-momentum states explored during the final doubling
#' iteration."
#' @param fn A function that returns the log of the posterior density. This
#' function should not return the negative log, following the notation of
#' Hoffman and Gelman (2014).
#' @param gr A function that returns a vector of gradients of the log of
#' the posterior density (same as with \code{fn}).
#' @param params.init A vector of initial parameter values.
#' @param Madapt The number of iterations during which to tune \code{eps}
#' if the dual averaging algorithm is used. Afterward the final \code{eps}
#' is used for the remaining iterations.
#' @param delta The target acceptance rate for the dual averaging
#' algorithm. See paper for interpretation of this for NUTS which does not
#' have a Metropolis step.
#' @param diagnostic Whether to return a list of diagnostic metrics about
#' the chain. Useful for assessing efficiency and tuning chain.
#' @param max_doublings Integer representing the maximum times the path
#' length should double within an MCMC iteration. Default of 8, so 256
#' steps.
#' @references Hoffman and Gelman (2014). The No-U-Turn sampler: Adaptively
#' setting path lengths in Hamiltonian Monte Carlo.
#' @return If \code{diagnostic} is FALSE (default), returns a matrix of
#' \code{nsim} samples from the posterior. Otherwise returns a list
#' containing samples ('par'),  vector of steps taken at each iteration
#' ('steps.taken'), and the total function and gradient
#' calls ('n.calls').
mcmc.nuts <- function(nsim, fn, gr, params.init, Madapt, eps=NULL,
                      delta=0.5, covar=NULL, diagnostic=FALSE, max_doublings=7){
    ## If using covariance matrix and Cholesky decomposition, redefine
    ## these functions to include this transformation. The algorithm will
    ## work in the transformed space
    if(!is.null(covar)){
        fn2 <- function(theta) fn(chd %*% theta)
        gr2 <- function(theta) gr(chd %*% theta)
        chd <- t(chol(covar))               # lower triangular Cholesky decomp.
        chd.inv <- solve(chd)               # inverse
        theta.cur <- chd.inv %*% params.init
    } else {
        fn2 <- fn; gr2 <- gr
        theta.cur <- params.init
    }
    theta.cur <- params.init
    theta.out <- matrix(NA, nrow=nsim, ncol=length(theta.cur))
    ## how many steps were taken at each iteration, useful for tuning
    j.results <- rep(NA, len=nsim)
    ## count the model calls as global variable; updated inside
    ## .buildtree. Some subtrees wont finish due to exit conditions so this
    ## is dynamic and not a simple formula like with HMC.
    assign("n.calls", value=0, envir=.GlobalEnv)
    useDA <- is.null(eps)               # whether to use DA algorithm
    if(useDA){
        ## Initialize the dual-averaging algorithm. Could make these arguments
        ## later.
        message(paste("No eps given so using dual averaging during first", Madapt, "steps."))
        epsvec <- Hbar <- epsbar <- rep(NA, length=Madapt+1)
        eps <- epsvec[1] <- find.epsilon(theta=theta.cur, fn=fn2, gr=gr2, eps=.1)
        mu <- log(10*eps)
        epsbar[1] <- 1; Hbar[1] <- 0; gamma <- 0.05; t0 <- 10; kappa <- 0.75
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
                                  fn=fn2, gr=gr2)
                theta.plus <- res$theta.plus
                r.plus <- res$r.plus
            } else {
                ## move in left direction
                res <- .buildtree(theta=theta.minus, r=r.minus, u=u, v=v,
                                  j=j, eps=eps, theta0=theta0, r0=r0,
                                  fn=fn2, gr=gr2)
                theta.minus <- res$theta.minus
                r.minus <- res$r.minus
            }
            ## test whether to accept this state
            if(is.na(res$s)) browser()
            if(res$s==1) {
                if(runif(n=1, min=0,max=1) <= res$n/n){
                    theta.cur <- res$theta.prime
                    theta.out[m,] <- res$theta.prime
                }
            }
            n <- n+res$n
            s <- res$s*.test.nuts(theta.plus, theta.minus, r.plus, r.minus)
            j <- j+1
            if(j>max_doublings & s) {warning("j larger than max_doublings, skipping to next m");break}
        }
        j.results[m] <- j-1
        if(useDA){
            ## Do the adapting of eps. Note that indexing is subtle here, the
            ## paper uses 0 but R needs to start at 1. I've thus offset
            ## every index by +1.
            if(m <= Madapt){
                Hbar[m+1] <-
                    (1-1/(m+t0))*Hbar[m] + (delta-res$alpha/res$nalpha)/(m+t0)
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
    message(paste0("NUTS diagnostics: Approximate leapfrog steps(min, median, max)=(",
                   paste(j.stats, collapse=","), ")"))
    message(paste("Total function calls:", n.calls))
    if(diagnostic){
        return(list(par=theta.out, steps.taken= 2^j.results,
                    n.calls=n.calls, epsvec=epsvec, epsbar=epsbar, Hbar=Hbar))
    } else {
        return(theta.out)
    }
}


#' Draw a slice sample for given position and momentum variables
.sample.u <- function(theta, r, fn)
    runif(n=1, min=0, max=exp(.calculate.H(theta=theta,r=r, fn=fn)))
#' Calculate the Hamiltonian value for position and momentum variables.
#'
#' @details This function currently assumes iid standard normal momentum
#' variables.
.calculate.H <- function(theta, r, fn) fn(theta)-(1/2)*sum(r^2)
#' Test whether a "U-turn" has occured in a branch of the binary tree
#' created by \ref\code{.buildtree} function.
.test.nuts <- function(theta.plus, theta.minus, r.plus, r.minus){
    theta.temp <- (theta.plus-theta.minus)
   as.numeric( theta.temp %*% r.minus >= 0 | theta.temp %*% r.plus >= 0)
}

#' A recursive function that builds a leapfrog trajectory using a balanced
#' binary tree.
#'
#' @references This is from the No-U-Turn sampler with dual averaging
#' (algorithm 6) of Hoffman and Gelman (2014).
#'
#' @details The function repeatedly doubles (in a random direction) until
#' either a U-turn occurs or the trajectory becomes unstable. This is the
#' 'efficient' version that samples uniformly from the path without storing
#' it. Thus the function returns a single proposed value and not the whole
#' trajectory.
#'
.buildtree <- function(theta, r, u, v, j, eps, theta0, r0, fn, gr,
                         delta.max=1000){
    if(j==0){
        ## base case, take one step in direction v
        eps <- v*eps
        r <- r+(eps/2)*gr(theta)
        theta <- theta+eps*r
        r <- r+(eps/2)*gr(theta)
        ## verify valid trajectory
        H <- .calculate.H(theta=theta, r=r, fn=fn)
        s <- H-log(u) + delta.max > 0
        n <- log(u) <= H
        ## ## Useful code for debugging. Returns entire path to global env.
        ## if(!exists('theta.trajectory'))
        ##     theta.trajectory <<- theta
        ## else
        ##     theta.trajectory <<- rbind(theta.trajectory, theta)
        temp <- .calculate.H(theta=theta, r=r, fn=fn)-
            .calculate.H(theta=theta0, r=r0, fn=fn)
        alpha <- min(exp(temp),1)
        if(exists('n.calls'))
            assign("n.calls", value=n.calls+5, envir=.GlobalEnv)
        else
            assign("n.calls", value=0, envir=.GlobalEnv)
        return(list(theta.minus=theta, theta.plus=theta, theta.prime=theta, r.minus=r,
                    r.plus=r, s=s, n=n, alpha=alpha, nalpha=1))
    } else {
        ## recursion - build left and right subtrees
        xx <- .buildtree(theta=theta, r=r, u=u, v=v, j=j-1, eps=eps,
                            theta0=theta0, r0=r0, fn=fn,gr=gr)
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
                                    fn=fn, gr=gr)
                theta.minus <- yy$theta.minus
                r.minus <- yy$r.minus
            } else {
                yy <- .buildtree(theta=theta.plus, r=r.plus, u=u, v=v,
                                    j=j-1, eps=eps, theta0=theta0, r0=r0,
                                    fn=fn, gr=gr)
                theta.plus <- yy$theta.plus
                r.plus <- yy$r.plus
            }
            ## This isn't in the paper but if both slice variables failed,
            ## then you get 0/0. So I skip this test
            nprime <- yy$n+ xx$n
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




#' Estimate a reasonable starting value for epsilon (step size) for a given
#' model, for use with Hamiltonian MCMC algorithms.
#'
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
find.epsilon <- function(theta,  fn, gr, eps=1){
    r <- rnorm(n=length(theta), mean=0, sd=1)
    ## Do one leapfrog step
    r.new <- r+(eps/2)*gr(theta)
    theta.new <- theta+eps*r.new
    r.new <- r.new+(eps/2)*gr(theta.new)
    H1 <- .calculate.H(theta=theta, r=r, fn=fn)
    H2 <- .calculate.H(theta=theta.new, r=r.new, fn=fn)
    a <- 2*(exp(H2)/exp(H1)>.5)-1
    k <- 1
    while( (exp(H2)/exp(H1))^a > 2^(-a) ){
        eps <- (2^a)*eps
        ## Do one leapfrog step
        r.new <- r+(eps/2)*gr(theta)
        theta.new <- theta+eps*r.new
        r.new <- r.new+(eps/2)*gr(theta.new)
        H2 <- .calculate.H(theta=theta.new, r=r.new, fn=fn)
        k <- k+1
        if(k>50) {
            warning("more than 50 iterations to find epsilon, stopping")
            break
        }
    }
    message(paste("Reasonable epsilon=", eps, "found after", k, "steps"))
    return(invisible(eps))
}
