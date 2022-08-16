## Copyright (C) 2013-2015 Kasper Kristensen
## License: GPL-2

##' Calculate one-step-ahead (OSA) residuals for a latent variable
##' model. (\emph{Beta version; may change without notice})
##'
##' Given a TMB latent variable model this function calculates OSA
##' standardized residuals that can be used for goodness-of-fit
##' assessment. The approach is based on a factorization of the joint
##' distribution of the \emph{observations} \eqn{X_1,...,X_n} into
##' successive conditional distributions.
##' Denote by
##' \deqn{F_n(x_n) = P(X_n \leq x_n | X_1 = x_1,...,X_{n-1}=x_{n-1} )}
##' the one-step-ahead CDF, and by
##' \deqn{p_n(x_n) = P(X_n = x_n | X_1 = x_1,...,X_{n-1}=x_{n-1} )}
##' the corresponding point probabilities (zero for continuous distributions).
##' In case of continuous observations the sequence
##' \deqn{\Phi^{-1}(F_1(X_1))\:,...,\:\Phi^{-1}(F_n(X_n))}
##' will be iid standard normal. These are referred to as the OSA residuals.
##' In case of discrete observations draw (unit) uniform variables
##' \eqn{U_1,...,U_n} and construct the randomized OSA residuals
##' \deqn{\Phi^{-1}(F_1(X_1)-U_1 p_1(X_1))\:,...,\:\Phi^{-1}(F_n(X_n)-U_n p_n(X_n))}
##' These are also iid standard normal.
##'
##' @section Choosing the method:
##' The user must specify the method used to calculate the residuals - see detailed list of method descriptions below.
##' We note that all the methods are based on approximations. While the default 'oneStepGaussianoffMode' often represents a good compromise between accuracy and speed, it cannot be assumed to work well for all model classes.
##' As a rule of thumb, if in doubt whether a method is accurate enough, you should always compare with the 'oneStepGeneric' which is considered the most accurate of the available methods.
##' \describe{
##' \item{method="fullGaussian"}{
##' This method assumes that the joint distribution of data \emph{and}
##' random effects is Gaussian (or well approximated by a
##' Gaussian). It does not require any changes to the user
##' template. However, if used in conjunction with \code{subset}
##' and/or \code{conditional} a \code{data.term.indicator} is required
##' - see the next method.
##' }
##' \item{method="oneStepGeneric"}{
##' This method calculates the one-step conditional probability
##' density as a ratio of Laplace approximations. The approximation is
##' integrated (and re-normalized for improved accuracy) using 1D
##' numerical quadrature to obtain the one-step CDF evaluated at each
##' data point. The method works in the continuous case as well as the
##' discrete case (\code{discrete=TRUE}).
##'
##' It requires a specification of a \code{data.term.indicator}
##' explained in the following. Suppose the template for the
##' observations given the random effects (\eqn{u}) looks like
##' \preformatted{
##'     DATA_VECTOR(x);
##'     ...
##'     nll -= dnorm(x(i), u(i), sd(i), true);
##'     ...
##' }
##'
##' Then this template can be augmented with a
##' \code{data.term.indicator = "keep"} by changing the template to
##' \preformatted{
##'     DATA_VECTOR(x);
##'     DATA_VECTOR_INDICATOR(keep, x);
##'     ...
##'     nll -= keep(i) * dnorm(x(i), u(i), sd(i), true);
##'     ...
##' }
##'
##' The new data vector (\code{keep}) need not be passed from \R. It
##' automatically becomes a copy of \code{x} filled with ones.
##'
##' Some extra parameters are essential for the method.
##' Pay special attention to the integration domain which must be specified either via \code{range} (continuous case) or \code{discreteSupport} (discrete case).
##' It may be useful to look at the one step predictive distributions on either log scale (\code{trace=2}) or natural scale (\code{trace=3}) to determine which alternative methods might be appropriate.
##' }
##' \item{method="oneStepGaussian"}{
##' This is a special case of the generic method where the one step
##' conditional distribution is approximated by a Gaussian (and can
##' therefore be handled more efficiently).
##' }
##' \item{method="oneStepGaussianOffMode"}{
##' This is an approximation of the "oneStepGaussian" method that
##' avoids locating the mode of the one-step conditional density.
##' }
##' \item{method="cdf"}{
##' The generic method can be slow due to the many function
##' evaluations used during the 1D integration (or summation in the
##' discrete case). The present method can speed up this process but
##' requires more changes to the user template. The above template
##' must be expanded with information about how to calculate the
##' negative log of the lower and upper CDF:
##' \preformatted{
##'     DATA_VECTOR(x);
##'     DATA_VECTOR_INDICATOR(keep, x);
##'     ...
##'     nll -= keep(i) * dnorm(x(i), u(i), sd(i), true);
##'     nll -= keep.cdf_lower(i) * log( pnorm(x(i), u(i), sd(i)) );
##'     nll -= keep.cdf_upper(i) * log( 1.0 - pnorm(x(i), u(i), sd(i)) );
##'     ...
##' }
##'
##' The specialized members \code{keep.cdf_lower} and
##' \code{keep.cdf_upper} automatically become copies of \code{x}
##' filled with zeros.
##' }
##' }
##'
##' @title Calculate one-step-ahead (OSA) residuals for a latent variable model.
##' @param obj Output from \code{MakeADFun}.
##' @param observation.name Character naming the observation in the template.
##' @param data.term.indicator Character naming an indicator data variable in the template (not required by all methods - see details).
##' @param method Method to calculate OSA (see details).
##' @param subset Index vector of observations that will be added one by one during OSA. By default \code{1:length(observations)} (with \code{conditional} subtracted).
##' @param conditional Index vector of observations that are fixed during OSA. By default the empty set.
##' @param discrete Are observations discrete? (assumed FALSE by default)
##' @param discreteSupport Possible outcomes of discrete distribution (\code{method="oneStepGeneric"} only).
##' @param range Possible range of the observations. (\code{method="oneStepGeneric"} only).
##' @param seed Randomization seed (discrete case only). If \code{NULL} the RNG seed is untouched by this routine (recommended for simulation studies).
##' @param parallel Run in parallel using the \code{parallel} package?
##' @param trace Logical; Trace progress? More options available for \code{method="oneStepGeneric"} - see details.
##' @param reverse Do calculations in opposite order to improve stability? (currently enabled by default for \code{oneStepGaussianOffMode} method only)
##' @param ... Control parameters for OSA method
##' @return \code{data.frame} with OSA \emph{standardized} residuals
##' in column \code{residual}. In addition, depending on the method, the output
##' includes selected characteristics of the predictive distribution (current row) given past observations (past rows), notably the \emph{conditional}
##' \describe{
##' \item{mean}{Expectation of the current observation}
##' \item{sd}{Standard deviation of the current observation}
##' \item{Fx}{CDF at current observation}
##' \item{px}{Density at current observation}
##' \item{nll}{Negative log density at current observation}
##' \item{nlcdf.lower}{Negative log of the lower CDF at current observation}
##' \item{nlcdf.upper}{Negative log of the upper CDF at current observation}
##' }
##' \emph{given past observations}.
##' @examples
##' ######################## Gaussian case
##' runExample("simple")
##' osa.simple <- oneStepPredict(obj, observation.name = "x", method="fullGaussian")
##' qqnorm(osa.simple$residual); abline(0,1)
##'
##' \dontrun{
##' ######################## Poisson case (First 100 observations)
##' runExample("ar1xar1")
##' osa.ar1xar1 <- oneStepPredict(obj, "N", "keep", method="cdf", discrete=TRUE, subset=1:100)
##' qqnorm(osa.ar1xar1$residual); abline(0,1)
##' }
oneStepPredict <- function(obj,
                           ## Names of data objects (not all are optional)
                           observation.name = NULL,
                           data.term.indicator = NULL,
                           method=c(
                               "oneStepGaussianOffMode",
                               "fullGaussian",
                               "oneStepGeneric",
                               "oneStepGaussian",
                               "cdf"),
                           subset = NULL,
                           conditional = NULL,
                           discrete = NULL,
                           discreteSupport = NULL,
                           range = c(-Inf, Inf),
                           seed = 123,
                           parallel = FALSE,
                           trace = TRUE,
                           reverse = (method == "oneStepGaussianOffMode"),
                           ...
                           ){
    if (missing(observation.name))
        stop("'observation.name' must define a data component")
    if (!(observation.name %in% names(obj$env$data)))
        stop("'observation.name' must be in data component")
    method <- match.arg(method)
    if (is.null(data.term.indicator)){
        if(method != "fullGaussian"){
            stop(paste0("method='",method,"' requires a 'data.term.indicator'"))
        }
    }
    if (!missing(discreteSupport) && !missing(range))
        stop("Cannot specify both 'discreteSupport' and 'range'")
    obs <- as.vector(obj$env$data[[observation.name]])
    if(is.null(discrete)){
        ndup <- sum(duplicated(obs))
        if(ndup > 0){
            warning("Observations do not look continuous. Number of duplicates = ", ndup)
            stop("Argument 'discrete' (TRUE/FALSE) must be specified.")
        }
        discrete <- FALSE
    } else {
        stopifnot(is.logical(discrete))
    }
    ## Using wrong method for discrete data ?
    if (discrete){
        if (! (method %in% c("oneStepGeneric", "cdf")) ){
            stop(paste0("method='",method,"' is not for discrete observations."))
        }
    }
    ## Default subset/permutation:
    if(is.null(subset)){
        subset <- 1:length(obs)
        subset <- setdiff(subset, conditional)
    }
    ## Check
    if(!is.null(conditional)){
        if(length(intersect(subset, conditional)) > 0){
            stop("'subset' and 'conditional' have non-empty intersection")
        }
    }
    unconditional <- setdiff(1:length(obs), union(subset, conditional))

    ## Args to construct copy of 'obj'
    args <- as.list(obj$env)[intersect(names(formals(MakeADFun)), ls(obj$env))]
    ## Use the best encountered parameter for new object
    if(length(obj$env$random))
        args$parameters <- obj$env$parList(par = obj$env$last.par.best)
    else
        args$parameters <- obj$env$parList(obj$env$last.par.best)
    ## Fix all non-random components of parameter list
    names.random <- unique(names(obj$env$par[obj$env$random]))
    names.all <- names(args$parameters)
    fix <- setdiff(names.all, names.random)
    map <- lapply(args$parameters[fix], function(x)factor(x*NA))
    ran.in.map <- names.random[names.random %in% names(args$map)]
    if(length(ran.in.map)) map <- c(map, args$map[ran.in.map]) # don't overwrite random effects mapping
    args$map <- map ## Overwrite map
    ## Find randomeffects character
    args$random <- names.random
    args$regexp <- FALSE
    ## Do we need to change, or take derivatives wrt., observations?
    ## (search the code to see if a method uses "observation(k,y)" or
    ## just "observation(k)").
    obs.are.variables <- (method != "cdf")
    ## Move data$name to parameter$name if necessary
    if (obs.are.variables) {
        args$parameters[observation.name] <- args$data[observation.name]
        args$data[observation.name] <- NULL
    }
    ## Make data.term.indicator in parameter list
    if(!is.null(data.term.indicator)){
        one <- rep(1, length(obs))
        zero <- rep(0, length(obs))
        if(method=="cdf"){
            args$parameters[[data.term.indicator]] <- cbind(one, zero, zero)
        } else {
            args$parameters[[data.term.indicator]] <- cbind(one)
        }
        ## Set attribute to tell the order of observations
        ord <- rep(NA, length(obs))
        ord[conditional] <- 0 ## First (out of bounds)
        ord[subset] <- seq_along(subset)
        ord[unconditional] <- length(obs) + 1 ## Never (out of bounds)
        if (any(is.na(ord))) {
            stop("Failed to determine the order of obervations")
        }
        attr(args$parameters[[data.term.indicator]], "ord") <- as.double(ord - 1)
    }
    ## Pretend these are *not observed*:
    if(length(unconditional)>0){
        if(is.null(data.term.indicator))
            stop("Failed to disable some data terms (because 'data.term.indicator' missing)")
        args$parameters[[data.term.indicator]][unconditional, 1] <- 0
    }
    ## Pretend these are *observed*:
    if(length(conditional)>0){
        if(is.null(data.term.indicator))
            stop("Failed to enable some data terms (because 'data.term.indicator' missing)")
        args$parameters[[data.term.indicator]][conditional, 1] <- 1
    }
    ## Make map for observations and indicator variables:
    makeFac <- function(x){
        fac <- as.matrix(x)
        fac[] <- 1:length(x)
        fac[conditional, ] <- NA
        fac[unconditional, ] <- NA
        fac[subset, ] <- 1:(length(subset)*ncol(fac)) ## Permutation
        factor(fac)
    }
    map <- list()
    if (obs.are.variables)
        map[[observation.name]] <- makeFac(obs)
    if(!is.null(data.term.indicator)){
        map[[data.term.indicator]] <- makeFac(args$parameters[[data.term.indicator]])
    }
    args$map <- c(args$map, map)
    ## New object be silent
    args$silent <- TRUE
    ## Create new object
    newobj <- do.call("MakeADFun", args)

    ## Helper function to loop through observations:
    nm <- names(newobj$par)
    obs.pointer <- which(nm == observation.name)

    if(method=="cdf"){
        tmp <- matrix(which(nm == data.term.indicator), ncol=3)
        data.term.pointer <- tmp[,1]
        lower.cdf.pointer <- tmp[,2]
        upper.cdf.pointer <- tmp[,3]
    } else {
        data.term.pointer <- which(nm == data.term.indicator)
        lower.cdf.pointer <- NULL
        upper.cdf.pointer <- NULL
    }
    observation <- local({
        obs.local <- newobj$par
        i <- 1:length(subset)
        function(k, y=NULL, lower.cdf=FALSE, upper.cdf=FALSE){
            ## Disable all observations later than k:
            obs.local[data.term.pointer[k<i]] <- 0
            ## On request, overwrite k'th observation:
            if(!is.null(y)) obs.local[obs.pointer[k]] <- y
            ## On request, get tail probs rather than point probs:
            if(lower.cdf | upper.cdf){
                obs.local[data.term.pointer[k]] <- 0
                if(lower.cdf) obs.local[lower.cdf.pointer[k]] <- 1
                if(upper.cdf) obs.local[upper.cdf.pointer[k]] <- 1
            }
            obs.local
        }
    })

    ## Parallel case: overload lapply
    if(parallel){
        ## mclapply uses fork => must set nthreads=1
        nthreads.restore <- TMB::openmp()
        on.exit( TMB::openmp( nthreads.restore ), add=TRUE)
        TMB::openmp(1)
        requireNamespace("parallel") # was library(parallel)
        lapply <- parallel::mclapply
    }
    ## Trace one-step functions
    tracefun <- function(k)if(trace)print(k)
    ## Apply a one-step method and generate common output assuming
    ## the method generates at least:
    ##   * nll
    ##   * nlcdf.lower
    ##   * nlcdf.upper
    applyMethod <- function(oneStepMethod){
        ord <- seq_along(subset)
        if (reverse) ord <- rev(ord)
        pred <- do.call("rbind", lapply(ord, oneStepMethod))
        pred <- as.data.frame(pred)[ord, ]
        pred$Fx <- 1 / ( 1 + exp(pred$nlcdf.lower - pred$nlcdf.upper) )
        pred$px <- 1 / ( exp(-pred$nlcdf.lower + pred$nll) +
                         exp(-pred$nlcdf.upper + pred$nll) )
        if(discrete){
            if(!is.null(seed)){
                ## Restore RNG on exit:
                Random.seed <- .GlobalEnv$.Random.seed
                on.exit(.GlobalEnv$.Random.seed <- Random.seed)
                set.seed(seed)
            }
            U <- runif(nrow(pred))
        } else {
            U <- 0
        }
        pred$residual <- qnorm(pred$Fx - U * pred$px)
        pred
    }

    ## ######################### CASE: oneStepGaussian
    if(method == "oneStepGaussian"){
        p <- newobj$par
        newobj$fn(p) ## Test eval
        oneStepGaussian <- function(k){
            tracefun(k)
            index <- subset[k]
            f <- function(y){
                newobj$fn(observation(k, y))
            }
            g <- function(y){
                newobj$gr(observation(k, y))[obs.pointer[k]]
            }
            opt <- nlminb(obs[index], f, g)
            H <- optimHess(opt$par, f, g)
            c(observation=obs[index], mean=opt$par, sd=sqrt(1/H))
        }
        ord <- seq_along(subset)
        if (reverse) ord <- rev(ord)
        pred <- do.call("rbind", lapply(ord, oneStepGaussian))
        pred <- as.data.frame(pred)[ord, ]
        pred$residual <- (pred$observation-pred$mean)/pred$sd
    }

    ## ######################### CASE: oneStepGaussianOffMode
    if(method == "oneStepGaussianOffMode"){
        p <- newobj$par
        newobj$fn(p) ## Test eval
        newobj$env$random.start <- expression({last.par[random]})
        oneStepGaussian <- function(k){
            tracefun(k)
            index <- subset[k]
            f <- function(y){
                newobj$fn(observation(k, y))
            }
            g <- function(y){
                newobj$gr(observation(k, y))[obs.pointer[k]]
            }
            c(observation=obs[index], nll = f(obs[index]), grad = g(obs[index]))
        }
        ord <- seq_along(subset)
        if (reverse) ord <- rev(ord)
        pred <- do.call("rbind", lapply(ord, oneStepGaussian))
        pred <- as.data.frame(pred)[ord, ]
        ################### Convert value and gradient to residual
        ## Need Lambert W function: x = W(x) * exp( W(x) ) , x > 0
        ## Vectorized in x and tested on extreme cases W(.Machine$double.xmin)
        ## and W(.Machine$double.xmax).
        W <- function(x){
            ## Newton: f(y)  = y * exp(y) - x
            ##         f'(y) = y * exp(y) + exp(y)
            rel.tol <- sqrt(.Machine$double.eps)
            logx <- log(x)
            fdivg <- function(y)(y - exp(logx - y)) / (1 + y)
            y <- pmax(logx, 0)
            while( any( abs( logx - log(y) - y) > rel.tol, na.rm=TRUE) ) {
                y <- y - fdivg(y)
            }
            y
        }
        getResid <- function(value, grad){
            Rabs <- sqrt( W( exp( 2*(value - log(sqrt(2*pi)) + log(abs(grad))) ) ) )
            R <- sign(grad) * Rabs
            R
        }
        nll0 <- newobj$fn(observation(0))
        R <- getResid( diff( c(nll0, pred$nll) ), pred$grad )
        M <- pred$observation - ifelse(pred$grad != 0, R * (R / pred$grad), 0)
        pred$mean <- M
        pred$residual <- R
    }

    ## ######################### CASE: oneStepGeneric
    if((method == "oneStepGeneric") && missing(discreteSupport)){
        p <- newobj$par
        newobj$fn(p) ## Test eval
        newobj$env$value.best <- -Inf ## <-- Never overwrite last.par.best
        nan2zero <- function(x)if(!is.finite(x)) 0 else x
        ## Set default configuration for this method (modify with '...'):
        formals(tmbprofile)$ytol <- 10  ## Tail tolerance (increase => more tail)
        formals(tmbprofile)$ystep <- .5 ## Grid spacing   (decrease => more accuracy)
        ## Handle discrete case
        if(discrete){
            formals(tmbprofile)$h <- 1
            integrate <- function(f, lower, upper, ...){
                grid <- ceiling(lower):floor(upper)
                list( value = sum( f(grid) ) )
            }
        }
        oneStepGeneric <- function(k){
            tracefun(k)
            ans <- try({
                index <- subset[k]
                f <- function(y){
                    newobj$fn(observation(k, y))
                }
                nll <- f(obs[index]) ## Marginal negative log-likelihood
                newobj$env$last.par.best <- newobj$env$last.par ## <-- used by tmbprofile
                slice <- tmbprofile(newobj, k, slice=TRUE,
                                    parm.range = range,...)
                spline <- splinefun(slice[[1]], slice[[2]])
                spline.range <- range(slice[[1]])
                if(trace >= 2){
                    plotfun <- function(slice, spline){
                        plot(slice, type="p", level=NULL)
                        plot(spline, spline.range[1], spline.range[2], add=TRUE)
                        abline(v=obs[index], lty="dashed")
                    }
                    if(trace >= 3){
                        slice$value <- exp( -(slice$value - nll) )
                        plotfun(slice, function(x)exp(-(spline(x) - nll)))
                    }
                    else
                        plotfun(slice, spline)
                }
                F1 <- integrate(function(x)exp(-(spline(x) - nll)),
                                spline.range[1],
                                obs[index])$value
                F2 <- integrate(function(x)exp(-(spline(x) - nll)),
                                obs[index] + discrete,
                                spline.range[2])$value
                mean <- integrate(function(x)exp(-(spline(x) - nll)) * x,
                                  spline.range[1],
                                  spline.range[2])$value / (F1 + F2)
                ## Was:
                ##  F1 <- integrate(Vectorize( function(x)nan2zero( exp(-(f(x) - nll)) ) ), -Inf, obs[index])$value
                ##  F2 <- integrate(Vectorize( function(x)nan2zero( exp(-(f(x) - nll)) ) ), obs[index], Inf)$value
                nlcdf.lower = nll - log(F1)
                nlcdf.upper = nll - log(F2)
                c(nll=nll, nlcdf.lower=nlcdf.lower, nlcdf.upper=nlcdf.upper, mean=mean)
            })
            if(is(ans, "try-error")) ans <- NaN
            ans
        }
        pred <- applyMethod(oneStepGeneric)
    }

    ## ######################### CASE: oneStepDiscrete
    if((method == "oneStepGeneric") && !missing(discreteSupport)){
        p <- newobj$par
        newobj$fn(p) ## Test eval
        obs <- as.integer(round(obs))
        if(is.null(discreteSupport)){
            warning("Setting 'discreteSupport' to ",min(obs),":",max(obs))
            discreteSupport <- min(obs):max(obs)
        }
        oneStepDiscrete <- function(k){
            tracefun(k)
            ans <- try({
                index <- subset[k]
                f <- function(y){
                    newobj$fn(observation(k, y))
                }
                nll <- f(obs[index]) ## Marginal negative log-likelihood
                F <- Vectorize(function(x)exp(-(f(x) - nll))) (discreteSupport)
                F1 <- sum( F[discreteSupport <= obs[index]] )
                F2 <- sum( F[discreteSupport >  obs[index]] )
                nlcdf.lower = nll - log(F1)
                nlcdf.upper = nll - log(F2)
                c(nll=nll, nlcdf.lower=nlcdf.lower, nlcdf.upper=nlcdf.upper)
            })
            if(is(ans, "try-error")) ans <- NaN
            ans
        }
        pred <- applyMethod(oneStepDiscrete)
    }

    ## ######################### CASE: fullGaussian
    if(method == "fullGaussian"){
        ## Same object with y random:
        args2 <- args
        args2$random <- c(args2$random, observation.name)
        ## Change map: Fix everything except observations
        fix <- data.term.indicator
        args2$map[fix] <- lapply(args2$map[fix],
                                 function(x)factor(NA*unclass(x)))
        newobj2 <- do.call("MakeADFun", args2)
        newobj2$fn() ## Test-eval to find mode
        mode <- newobj2$env$last.par
        GMRFmarginal <- function (Q, i, ...) {
            ind <- 1:nrow(Q)
            i1 <- (ind)[i]
            i0 <- setdiff(ind, i1)
            if (length(i0) == 0)
                return(Q)
            Q0 <- as(Q[i0, i0, drop = FALSE], "symmetricMatrix")
            L0 <- Cholesky(Q0, ...)
            ans <- Q[i1, i1, drop = FALSE] - Q[i1, i0, drop = FALSE] %*%
                solve(Q0, Q[i0, i1, drop = FALSE])
            ans
        }
        h <- newobj2$env$spHess(mode, random=TRUE)
        i <- which(names(newobj2$env$par[newobj2$env$random]) == observation.name)
        Sigma <- solve( as.matrix( GMRFmarginal(h, i) ) )
        res <- obs[subset] - mode[i]
        L <- t(chol(Sigma))
        pred <- data.frame(residual = as.vector(solve(L, res)))
    }

    ## ######################### CASE: cdf
    if(method == "cdf"){
        p <- newobj$par
        newobj$fn(p) ## Test eval
        newobj$env$random.start <- expression(last.par[random])
        cdf <- function(k){
            tracefun(k)
            nll <- newobj$fn(observation(k))
            lp <- newobj$env$last.par
            nlcdf.lower <- newobj$fn(observation(k, lower.cdf = TRUE))
            newobj$env$last.par <- lp ## restore
            nlcdf.upper <- newobj$fn(observation(k, upper.cdf = TRUE))
            newobj$env$last.par <- lp ## restore
            c(nll=nll, nlcdf.lower=nlcdf.lower, nlcdf.upper=nlcdf.upper)
        }
        pred <- applyMethod(cdf)
    }

    pred
}

## Goodness of fit residuals based on an approximate posterior
## sample. (\emph{Beta version; may change without notice})
##
## Denote by \eqn{(u, x)} the pair of the true un-observed random effect
## and the data. Let a model specification be given in terms of the
## estimated parameter vector \eqn{\theta} and let \eqn{u^*} be a
## sample from the conditional distribution of \eqn{u} given
## \eqn{x}. If the model specification is correct, it follows that the
## distribution of the pair \eqn{(u^*, x)} is the same as the distribution
## of \eqn{(u, x)}. Goodness-of-fit can thus be assessed by proceeding as
## if the random effect vector were observed, i.e check that \eqn{u^*}
## is consistent with prior model of the random effect and that \eqn{x}
## given \eqn{u^*} agrees with the observation model.
##
## This function can carry out the above procedure for many TMB models
## under the assumption that the true posterior is well approximated by a
## Gaussian distribution.
##
## First a draw from the Gaussian posterior distribution \eqn{u^*} is
## obtained based on the mode and Hessian of the random effects given the
## data.
## This sample uses sparsity of the Hessian and will thus work for large systems.
##
## An automatic standardization of the sample can be carried out \emph{if
## the observation model is Gaussian} (\code{fullGaussian=TRUE}). In this
## case the prior model is obtained by disabling the data term and
## calculating mode and Hessian. A \code{data.term.indicator} must be
## given in order for this to work. Standardization is performed using
## the sparse Cholesky of the prior precision.
## By default, this step does not use a fill reduction permutation \code{perm=FALSE}.
## This is often superior wrt. to interpretation of the.
## the natural order of the parameter vector is used \code{perm=FALSE}
## which may be superior wrt. to interpretation. Otherwise
## \code{perm=TRUE} a fill-reducing permutation is used while
## standardizing.
## @references Waagepetersen, R. (2006). A Simulation-based Goodness-of-fit Test for Random Effects in Generalized Linear Mixed Models. Scandinavian journal of statistics, 33(4), 721-731.
## @param obj TMB model object from \code{MakeADFun}.
## @param observation.name Character naming the observation in the template.
## @param data.term.indicator Character naming an indicator data variable in the template. Only used if \code{standardize=TRUE}.
## @param standardize Logical; Standardize sample with the prior covariance ? Assumes all latent variables are Gaussian.
## @param as.list Output posterior sample, and the corresponding standardized residual, as a parameter list ?
## @param perm Logical; Use a fill-reducing ordering when standardizing ?
## @param fullGaussian Logical; Flag to signify that the joint distribution of random effects and data is Gaussian.
## @return List with components \code{sample} and \code{residual}.
oneSamplePosterior <- function(obj,
                               observation.name = NULL,
                               data.term.indicator = NULL,
                               standardize = TRUE,
                               as.list = TRUE,
                               perm = FALSE,
                               fullGaussian = FALSE){
    ## Draw Gaussian posterior sample
    tmp <- obj$env$MC(n=1, keep=TRUE, antithetic=FALSE)
    samp <- as.vector( attr(tmp, "samples") )
    ## If standardize
    resid <- NULL
    if (standardize) {
        ## Args to construct copy of 'obj'
        args <- as.list(obj$env)[intersect(names(formals(MakeADFun)), ls(obj$env))]
        ## Use the best encountered parameter for new object
        args$parameters <- obj$env$parList(par = obj$env$last.par.best)
        ## Make data.term.indicator in parameter list
        obs <- obj$env$data[[observation.name]]
        nobs <- length(obs)
        zero <- rep(0, nobs)
        if ( ! fullGaussian )
            args$parameters[[data.term.indicator]] <- zero
        ## Fix all non-random components of parameter list
        names.random <- unique(names(obj$env$par[obj$env$random]))
        names.all <- names(args$parameters)
        fix <- setdiff(names.all, names.random)
        map <- lapply(args$parameters[fix], function(x)factor(x*NA))
        args$map <- map ## Overwrite map
        ## If 'fullGaussian == TRUE' turn 'obs' into a random effect
        if (fullGaussian) {
            names.random <- c(names.random, observation.name)
            args$parameters[[observation.name]] <- obs
        }
        ## Find randomeffects character
        args$random <- names.random
        args$regexp <- FALSE
        ## New object be silent
        args$silent <- TRUE
        ## Create new object
        newobj <- do.call("MakeADFun", args)
        ## Construct Hessian and Cholesky
        newobj$fn()
        ## Get Cholesky and prior mean
        ## FIXME: We are using the mode as mean. Consider skewness
        ## correction similar to 'bias.correct' in 'sdreport'.
        L <- newobj$env$L.created.by.newton
        mu <- newobj$env$last.par
        ## If perm == FALSE redo Cholesky with natural ordering
        if ( ! perm ) {
            Q <- newobj$env$spHess(mu, random=TRUE)
            L <- Matrix::Cholesky(Q, super=TRUE, perm=FALSE)
        }
        ## If 'fullGaussian == TRUE' add 'obs' to the sample
        if (fullGaussian) {
            tmp <- newobj$env$par * NA
            tmp[names(tmp) == observation.name] <- obs
            tmp[names(tmp) != observation.name] <- samp
            samp <- tmp
        }
        ## Standardize ( P * Q * P^T = L * L^T )
        r <- samp - mu
        rp <- r[L@perm + 1]
        Lt <- Matrix::t(
            as(L, "sparseMatrix")
        )
        resid <- as.vector( Lt %*% rp )
    }
    if (as.list) {
        if (standardize) obj <- newobj
        par <- obj$env$last.par.best
        asList <- function(samp) {
            par[obj$env$random] <- samp
            samp <- obj$env$parList(par=par)
            nm <- unique(names(obj$env$par[obj$env$random]))
            samp[nm]
        }
        samp <- asList(samp)
        if (!is.null(resid))
            resid <- asList(resid)
    }
    ans <- list()
    ans$sample <- samp
    ans$residual <- resid
    ans
}

if(FALSE) {
    library(TMB)
    runExample("MVRandomWalkValidation", exfolder="../../tmb_examples/validation")
    set.seed(1)
    system.time( qw <- TMB:::oneSamplePosterior(obj, "obs", "keep") )
    qqnorm(as.vector(qw$residual$u)); abline(0,1)
    runExample("rickervalidation", exfolder="../../tmb_examples/validation")
    set.seed(1)
    system.time( qw <- TMB:::oneSamplePosterior(obj, "Y", "keep") )
    qqnorm(as.vector(qw$residual$X)); abline(0,1)
    runExample("ar1xar1")
    set.seed(1)
    system.time( qw <- TMB:::oneSamplePosterior(obj, "N", "keep") )
    qqnorm(as.vector(qw$residual$eta)); abline(0,1)
}
