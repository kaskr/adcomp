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
#' @param ... Further arguments to be passed to the algorithm. See help
#' files for the samplers for further arguments.
#' @example inst/examples/mcmc_examples.R
mcmc <- function(obj, nsim, algorithm, params.init=NULL, ...){
    ## Initialization for all algorithms
    algorithm <- match.arg(algorithm, choices=c("HMC", "NUTS"))
    fn <- function(x) {
        z <- -obj$fn(x)
        if(is.nan(z)){
            warning(paste("replacing NaN w/ Inf at:", paste(x, collapse=" ")))
            z <- Inf
        }
        return(z)
    }
    gr <- function(x) -as.vector(obj$gr(x))
    obj$env$beSilent()                  # silence console output
    ## argument checking
    if(is.null(params.init)){
        params.init <- obj$par
    } else if(length(params.init) != length(obj$par)){
        stop("params.init is wrong length")
    }
    ## Select and run the chain.
    if(algorithm=="HMC")
        mcmc.out <-
            mcmc.hmc(nsim=nsim, fn=fn, gr=gr, params.init=params.init, ...)
    else if(algorithm=="NUTS")
        mcmc.out <-
            mcmc.nuts(nsim=nsim, fn=fn, gr=gr, params.init=params.init, ...)
    ## Clean up returned matrix
    mcmc.out <- as.data.frame(mcmc.out)
    names(mcmc.out) <- names(obj$par)
    return(invisible(mcmc.out))
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
#' @references Neal, R. M. 2011. MCMC using Hamiltonian dynamics.
#' @return A matrix of \code{nsim} samples from the posterior.
mcmc.hmc <- function(nsim, L, eps, fn, gr, params.init){
    theta.cur <- params.init
    accepted <- rep(NA, nsim)
    theta.out <- matrix(NA, nrow=nsim, ncol=length(params.init))
    for(m in 1:nsim){
        r.cur <- r.new <- rnorm(length(params.init),0,1)
        theta.new <- theta.cur
        theta.leapfrog <- matrix(NA, nrow=L, ncol=length(theta.cur))
        r.leapfrog <- matrix(NA, nrow=L, ncol=length(theta.cur))
        ## Make a half step for first iteration
        r.new <- r.new+eps*gr(theta.new)/2
        for(i in 1:L){
            theta.leapfrog[i,] <- theta.new
            r.leapfrog[i,] <- r.new
            theta.new <- theta.new+eps*r.new
            ## Full step except at end
            if(i!=L) r.new <- r.new+eps*gr(theta.new)
        }
        ## half step for momentum at the end
        r.new <- r.new+eps*gr(theta.new)/2
        ## negate r to make proposal symmetric
        r.new <- -r.new
        if(runif(1) <
           exp(-fn(theta.cur)+fn(theta.new)+ sum(r.cur^2)/2-sum(r.new^2)/2)){
            ## accept the proposed state
            theta.cur <- theta.new
            accepted[m] <- TRUE
        } else {
            ## otherwise reject it and stay there
            accepted[m] <- FALSE
        }
        theta.out[m,] <- theta.cur
    }
    message(paste("Acceptance rate = ", round(mean(accepted),1)))
    return(theta.out)
}

#' [BETA VERSION] Sample from a posterior using the No-U-Turn sampler.
#'
#' @references This is from 'efficient' No-U-Turn sampler (algorithm 3) of
#' Hoffman and Gelman (2014).
.buildtree <- function(theta, r, u, v, j, eps, fn, gr, delta.max=1000){
    if(j==0){
        ## base case, take one step in direction v
        eps <- v*eps
        r <- r+(eps/2)*gr(theta)
        theta <- theta+eps*r
        r <- r+(eps/2)*gr(theta)
        ## verify valid trajectory
        H <- .calculate.H(theta=theta, r=r, fn=fn)
        s <- H-log(u) + delta.max > 0
        slice <- log(u) <= H
        ## if(!s) print(paste("invalid s at k=", k))
        return(list(theta.minus=theta, theta.plus=theta, theta.prime=theta, r.minus=r,
                    r.plus=r, s=s, slice=slice))
    } else {
        ## recursion - build left and right subtrees
        xx <- .buildtree(theta=theta, r=r, u=u, v=v, j=j-1, eps=eps,
                                  fn=fn,gr=gr)
        theta.minus <- xx$theta.minus
        theta.plus <- xx$theta.plus
        theta.prime <- xx$theta.prime
        r.minus <- xx$r.minus
        r.plus <- xx$r.plus
        slice.new <- xx$slice
        s <- xx$s
        if(xx$s==1){
            if(v== -1){
                yy <-
                    .buildtree(theta=theta.minus, r=r.minus, u=u, v=v,
                                        j=j-1, eps=eps, fn=fn, gr=gr)
                theta.minus <- yy$theta.minus
                r.minus <- yy$r.minus
            } else {
                yy <-
                    .buildtree(theta=theta.plus, r=r.plus, u=u, v=v,
                                        j=j-1, eps=eps, fn=fn, gr=gr)
                theta.plus <- yy$theta.plus
                r.plus <- yy$r.plus
            }
            ## This isn't in the paper but if both slice variables failed,
            ## then you get 0/0. So I skip this test
            temp <- yy$slice+ xx$slice
            if(temp!=0){
            ## choose whether to keep this theta
                if(runif(n=1, min=0, max=1) <= yy$slice/temp)
                    theta.prime <- yy$theta.prime
            }
            ## print(yy$slice/(yy$slice+ xx$slice))
            slice.new <- xx$slice + yy$slice
            ## check for valid proposal
            test <- .test.nuts(theta.plus=theta.plus,
                              theta.minus=theta.minus, r.plus=r.plus,
                              r.minus=r.minus)
            ## if(!test) print(paste("U turn at k=", k))
            ## check if any of the stopping conditions were met
            s <- xx$s*yy$s*test
        }
        return(list(theta.minus=theta.minus, theta.plus=theta.plus,
                    theta.prime=theta.prime,
                    r.minus=r.minus, r.plus=r.plus, s=s, slice=slice.new))
    }
}

#' Draw MCMC samples from a model posterior using a Hamiltonian sampler.
#'
#' @param nsim The number of samples to return.
#' @param eps The length of the leapfrog steps.
#' @param fn A function that returns the log of the posterior density. This
#' function should not return the negative log, following the notation of
#' Hoffman and Gelman (2014).
#' @param gr A function that returns a vector of gradients of the log of
#' the posterior density (same as with \code{fn}).
#' @param params.init A vector of initial parameter values.
#' @references Neal, R. M. 2011. MCMC using Hamiltonian dynamics.
#' @return A matrix of \code{nsim} samples from the posterior.

mcmc.nuts <- function(nsim, fn, gr, params.init, eps){
    theta.cur <- params.init
    theta.out <- matrix(NA, nrow=nsim, ncol=length(params.init))
    j.results <- rep(NA, len=nsim)      # For diagnostics
    for(m in 1:nsim){
        ## initialize
        theta.out[m,] <- theta.minus <- theta.plus <- theta.cur
        r.cur <- r.plus <- r.minus <- rnorm(length(theta.cur),0,1)
        ## Draw a slice variable u
        u <- .sample.u(theta=theta.cur, r=r.cur, fn=fn)
        j <- 0; n <- 1; s <- 1
        while(s==1) {
            v <- sample(x=c(1,-1), size=1)
            if(v==1){
                ## move in right direction
                res <- .buildtree(theta=theta.plus, r=r.plus, u=u, v=v,
                                       j=j, eps=eps, fn=fn, gr=gr)
                theta.plus <- res$theta.plus
                r.plus <- res$r.plus
            } else {
                ## move in left direction
                res <- .buildtree(theta=theta.minus, r=r.minus, u=u, v=v,
                                       j=j, eps=eps, fn=fn, gr=gr)
                theta.minus <- res$theta.minus
                r.minus <- res$r.minus
            }
            ## test whether to accept this state
            if(res$s==1) {
                if(runif(n=1, min=0,max=1) <= res$n/n){
                    theta.cur <- res$theta.prime
                    theta.out[m,] <- res$theta.prime
                }
            }
            n <- n+res$n
            s <- res$s*.test.nuts(theta.plus, theta.minus, r.plus, r.minus)
            j <- j+1
            if(j>10 & s) {warning("j got to >10, skipping to next m");break}
        }
        j.results[m] <- j-1
    }
    j.stats <- 2^(c(min(j.results), median(j.results), max(j.results)))
    message(paste0("NUTS diagnostics: Approximate leapfrog steps(min, median, max)=(",
                   paste(j.stats, collapse=","), ")"))
    return(theta.out)
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
#' @references This is from 'efficient' No-U-Turn sampler (algorithm 3) of
#' Hoffman and Gelman (2014).
#'
#' @details The function repeatedly doubles (in a random direction) until
#' either a U-turn occurs or the trajectory becomes unstable.
#'
.buildtree <- function(theta, r, u, v, j, eps, fn, gr, delta.max=1000){
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
        ## if(!s) print(paste("invalid s at k=", k))
        return(list(theta.minus=theta, theta.plus=theta, theta.prime=theta, r.minus=r,
                    r.plus=r, s=s, n=n))
    } else {
        ## recursion - build left and right subtrees
        xx <- .buildtree(theta=theta, r=r, u=u, v=v, j=j-1, eps=eps,
                                  fn=fn,gr=gr)
        theta.minus <- xx$theta.minus
        theta.plus <- xx$theta.plus
        theta.prime <- xx$theta.prime
        r.minus <- xx$r.minus
        r.plus <- xx$r.plus
        s <- xx$s
        n.new <- xx$n
        ## If it didn't fail, update the above quantities
        if(s==1){
            if(v== -1){
                yy <-
                    .buildtree(theta=theta.minus, r=r.minus, u=u, v=v,
                                        j=j-1, eps=eps, fn=fn, gr=gr)
                theta.minus <- yy$theta.minus
                r.minus <- yy$r.minus
            } else {
                yy <-
                    .buildtree(theta=theta.plus, r=r.plus, u=u, v=v,
                                        j=j-1, eps=eps, fn=fn, gr=gr)
                theta.plus <- yy$theta.plus
                r.plus <- yy$r.plus
            }
            ## This isn't in the paper but if both slice variables failed,
            ## then you get 0/0. So I skip this test
            temp <- yy$n+ xx$n
            if(temp!=0){
            ## choose whether to keep this theta
                if(runif(n=1, min=0, max=1) <= yy$n/temp)
                    theta.prime <- yy$theta.prime
            }
            n.new <- xx$n + yy$n
            ## check for valid proposal
            test <- .test.nuts(theta.plus=theta.plus,
                              theta.minus=theta.minus, r.plus=r.plus,
                              r.minus=r.minus)
            ## if(!test) print(paste("U turn at k=", k))
            ## check if any of the stopping conditions were met
            s <- xx$s*yy$s*test
        }
        return(list(theta.minus=theta.minus, theta.plus=theta.plus,
                    theta.prime=theta.prime,
                    r.minus=r.minus, r.plus=r.plus, s=s, n=n.new))
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





