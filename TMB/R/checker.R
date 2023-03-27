##' Check consistency of various parts of a TMB implementation.
##' Requires that user has implemented simulation code for the data and
##' optionally random effects. (\emph{Beta version; may change without
##' notice})
##'
##' This function checks that the simulation code of random effects and
##' data is consistent with the implemented negative log-likelihood
##' function. It also checks whether the approximate \emph{marginal}
##' score function is central indicating whether the Laplace
##' approximation is suitable for parameter estimation.
##'
##' Denote by \eqn{u} the random effects, \eqn{\theta} the parameters
##' and by \eqn{x} the data.  The main assumption is that the user has
##' implemented the joint negative log likelihood \eqn{f_{\theta}(u,x)}
##' satisfying
##' \deqn{\int \int \exp( -f_{\theta}(u,x) ) \:du\:dx = 1}
##' It follows that the joint and marginal score functions are central:
##' \enumerate{
##'   \item \eqn{E_{u,x}\left[\nabla_{\theta}f_{\theta}(u,x)\right]=0}
##'   \item \eqn{E_{x}\left[\nabla_{\theta}-\log\left( \int \exp(-f_{\theta}(u,x))\:du \right) \right]=0}
##' }
##' For each replicate of \eqn{u} and \eqn{x} joint and marginal
##' gradients are calculated. Appropriate centrality tests are carried
##' out by \code{\link{summary.checkConsistency}}.  An asymptotic
##' \eqn{\chi^2} test is used to verify the first identity. Power of
##' this test increases with the number of simulations \code{n}.  The
##' second identity holds \emph{approximately} when replacing the
##' marginal likelihood with its Laplace approximation. A formal test
##' would thus fail eventually for large \code{n}. Rather, the gradient
##' bias is transformed to parameter scale (using the estimated
##' information matrix) to provide an estimate of parameter bias caused
##' by the Laplace approximation.
##'
##' @section Simulation/re-estimation:
##' A full simulation/re-estimation study is performed when \code{estimate=TRUE}.
##' By default \link[stats]{nlminb} will be used to perform the minimization, and output is stored in a separate list component 'estimate' for each replicate.
##' Should a custom optimizer be needed, it can be passed as a user function via the same argument (\code{estimate}).
##' The function (\code{estimate}) will be called for each simulation as \code{estimate(obj)} where \code{obj} is the simulated model object.
##' Current default corresponds to \code{estimate = function(obj) nlminb(obj$par,obj$fn,obj$gr)}.
##' @title Check consistency and Laplace accuracy
##' @param obj Object from \code{MakeADFun}
##' @param par Parameter vector (\eqn{\theta}) for simulation. If
##'     unspecified use the best encountered parameter of the object.
##' @param hessian Calculate the hessian matrix for each replicate ?
##' @param estimate Estimate parameters for each replicate ?
##' @param n Number of simulations
##' @param observation.name Optional; Name of simulated observation
##' @return List with gradient simulations (joint and marginal)
##' @seealso \code{\link{summary.checkConsistency}}, \code{\link{print.checkConsistency}}
##' @examples
##' \dontrun{
##' runExample("simple")
##' chk <- checkConsistency(obj)
##' chk
##' ## Get more details
##' s <- summary(chk)
##' s$marginal$p.value  ## Laplace exact for Gaussian models }
checkConsistency <- function(obj,
                             par = NULL,
                             hessian = FALSE,
                             estimate = FALSE,
                             n = 100,
                             observation.name = NULL
                             ) {
    ## Optimizer
    if (!is.logical(estimate)) {
        Optimizer <- match.fun(estimate)
        estimate <- TRUE
    } else {
        Optimizer <- function(obj) nlminb(obj$par, obj$fn, obj$gr)
    }
    ## Args to construct copy of 'obj'
    args <- as.list(obj$env)[intersect(names(formals(MakeADFun)), ls(obj$env))]
    ## Determine parameter and full parameter to use
    r0 <- r <- obj$env$random
    if( is.null(par) ) {
        ## Default case: Optimization has been carried out by user
        if (is.null(obj$env$last.par.best)) {
            stop("'par' not specified.")
        }
        parfull <- obj$env$last.par.best
        if( any(r) ) par <- parfull[-r] else par <- parfull
    } else {
        ## Custom case: User specifies parameter vector (fixed effects)
        parfull <- obj$env$par
        if( any(r) ) parfull[-r] <- par else parfull <- par
    }
    ## Get names of random effects (excluding profiled parameters)
    if(any(obj$env$profile)) {
        r0 <- r[ ! as.logical(obj$env$profile) ]
        names.profile <- unique(names(parfull[r[as.logical(obj$env$profile)]]))
    } else {
        names.profile <- NULL
    }
    names.random <- unique(names(parfull[r0]))
    ## Use 'parfull' for new object
    args$parameters <- obj$env$parList(par, par = parfull)
    ## Fix all profiled parameters
    map.profile <- lapply(args$parameters[names.profile], function(x)factor(x*NA))
    args$map <- c(args$map, map.profile)
    ## Find randomeffects character
    args$random <- names.random
    args$regexp <- FALSE
    ## Are we in 'fast' (no retape) mode ?
    fast <- !is.null(observation.name)
    if (fast) {
        ## Move data -> parameters
        ## Note: We really do need to know 'observation.name'. There
        ## could be other (deterministic) items in 'data'...
        args$parameters <- c(args$data[observation.name], args$parameters)
        args$data[observation.name] <- NULL
    }
    ## Create new object
    newobj <- do.call("MakeADFun", args)
    newobj0 <- newobj ## backup
    if (fast) {
        parobs <- names(newobj$par) %in% observation.name
        ## NOTE: Simulation is stored as part of 'newobj$env$par'
        expandpar <- function(par) {
            ## Incudes par fixed *and* simulation:
            ans <- newobj0$env$par[newobj0$env$lfixed()]
            ans[!parobs] <- par
            ans
        }
        ## FIXME: No 'obj$he()' in this object
        newobj <- list(fn=function(x)newobj0$fn(expandpar(x)),
                       gr=function(x)newobj0$gr(expandpar(x))[!parobs],
                       par=newobj0$par[!parobs],
                       env=newobj0$env
                       )
    }
    doSim <- function(...) {
        simdata <- newobj0$simulate(newobj0$env$par, complete=TRUE)
        if (!fast) {
            newobj$env$data <- simdata
        }
        ## Check that random effects have been simulated
        haveRandomSim <- all( names.random %in% names(simdata) )
        ## Set good inner starting values
        if (haveRandomSim) {
            if (fast) {
                for (nm in names.random) {
                    newobj$env$par[names(newobj$env$par) == nm] <- simdata[[nm]]
                }
            } else {
                newobj$env$parameters[names.random] <- simdata[names.random]
            }
        }
        ## FIXME: Mapped random effects not supported (yet) for 'fast' approach
        if (!fast && haveRandomSim) {
            ## Snippet taken from MakeADFun to account for mapped parameters:
            map <- args$map[names(args$map) %in% names.random]
            if (length(map) > 0) {
                param.map <- lapply(names(map), function(nam) {
                    updateMap(newobj$env$parameters[[nam]], map[[nam]])
                })
                keepAttrib(newobj$env$parameters[names(map)]) <- param.map
            }
        }
        if (fast) {
            ## Set simulated data
            for (nm in observation.name) {
                newobj$env$par[names(newobj$env$par) == nm] <- simdata[[nm]]
            }
            ## Set inits
            newobj$env$last.par.best <- newobj$env$par
            newobj$env$value.best <- Inf
        } else {
            ## This approach *must* redo Cholesky
            newobj$env$L.created.by.newton <- NULL
            newobj$env$retape()
        }
        ans <- list()
        if (haveRandomSim) {
            ans$gradientJoint <- newobj$env$f(order=1)
            if(!is.null(newobj$env$random))
                ans$gradientJoint <- ans$gradientJoint[-newobj$env$random]
            if (fast)
                ans$gradientJoint <- ans$gradientJoint[!parobs]
        }
        ans$gradient <- newobj$gr(par)
        if (hessian) ans$hessian <- optimHess(par, newobj$fn, newobj$gr)
        if (estimate) {
            newobj$par <- par ## just in case...
            ans$objective.true <- newobj$fn(par)
            ans$estimate <- try(Optimizer(newobj))
        }
        ans
    }
    ans <- lapply(seq_len(n), doSim)
    attr(ans, "par") <- par
    class(ans) <- "checkConsistency"
    ans
}

##' Summarize output from \code{\link{checkConsistency}}
##'
##' @title Summarize output from \code{\link{checkConsistency}}
##' @param object Output from \code{\link{checkConsistency}}
##' @param na.rm Logical; Remove failed simulations ?
##' @param ... Not used
##' @return List of diagnostics
##' @method summary checkConsistency
##' @S3method summary checkConsistency
summary.checkConsistency <- function(object, na.rm=FALSE, ...) {
    ans <- list()
    ans$par <- attr(object, "par")
    getMat <- function(name) {
        do.call("cbind",
                lapply(object,
                       function(x)
                           as.vector(x[[name]])))
    }
    ans$gradientJoint <- getMat( "gradientJoint" )
    ans$gradient      <- getMat( "gradient" )
    ## Check simulation
    check <- function(mat) {
        if(!is.matrix(mat)) return( list(p.value=NA, bias=NA) )
        if (na.rm) {
            fail <- as.logical( colSums( !is.finite(mat) ) )
            mat <- mat[, !fail, drop=FALSE]
        }
        mu <- rowMeans(mat)
        npar <- length(mu)
        nsim <- ncol(mat)
        bias <- p.value <- NULL
        if(nsim < npar) {
            stop("Too few simulations ", nsim, " compared to number of parameters ", npar)
        }
        ## Variance of score = Information
        H <- var(t(mat))
        iH <- try(solve(H), silent=TRUE)
        if(is(iH, "try-error")) {
            warning("Failed to invert information matrix")
            bias <- attr(object, "par") * NA
            p.value <- NA
        } else {
            mu.scaled <- sqrt(nsim) * mu
            q <- as.vector( t(mu.scaled) %*% iH %*% mu.scaled )
            p.value <- 1 - pchisq(q, df=npar)
            bias <- -iH %*% mu
        }
        bias <- as.vector(bias)
        names(bias) <- names(attr(object, "par"))
        list(p.value=p.value, bias=bias)
    }
    ans$joint <- check( ans$gradientJoint )
    ans$marginal <- check( ans$gradient )
    ## Simulation study
    have.estimate <- !is.null(object[[1]]$estimate)
    if (have.estimate) {
        getEstMat <- function(name) {
            do.call("cbind",
                    lapply(object,
                           function(x)
                               as.vector(x$estimate[[name]])))
        }
        est <- list()
        est$par <- t(getEstMat("par"))
        colnames(est$par) <- names(ans$par)
        est$par <- as.data.frame(est$par)
        est$objective <- drop(getEstMat("objective"))
        est$deviance <- 2 * ( drop(getMat("objective.true")) - est$objective )
        est$deviance.p.value <-
            ks.test(est$deviance, "pchisq", df = length(ans$par))$p.value
        ans$convergence <- drop(getEstMat("convergence"))
        ## Set it
        ans$estimate <- est
    }
    ans
}

##' Print diagnostics output from \code{\link{checkConsistency}}
##'
##' @title Print output from \code{\link{checkConsistency}}
##' @param x Output from \code{\link{checkConsistency}}
##' @param ... Not used
##' @return NULL
##' @method print checkConsistency
##' @S3method print checkConsistency
print.checkConsistency <- function(x, ...) {
    s <- summary(x)
    cat("Parameters used for simulation:\n")
    print(s$par)
    cat("\n")
    cat("Test correct simulation (p.value):\n")
    print(s$joint$p.value)
    alpha <- .05 ## FIXME: Perhaps make option
    s$sim.ok <- ( s$joint$p.value > alpha )
    if(is.na(s$sim.ok))
        cat("Full simulation was not available\n")
    else if(!s$sim.ok)
        cat("Simulation does *not* appear to be correct !!!\n")
    else
        cat("Simulation appears to be correct\n")
    ## Check Laplace:
    cat("\n")
    cat("Estimated parameter bias:\n")
    print(s$marginal$bias)
    ## Estimate info:
    if (!is.null(s$estimate)) {
        cat("\n")
        cat("summary(.)$estimate contains:\n")
        print(names(s$estimate))
    }
    invisible(x)
}

if(FALSE) {
    library(TMB)
    runExample("sam", exfolder="../../tmb_examples")
    set.seed(123)
    qw <- checkConsistency(obj, opt$par, n=100)
    print.checkConsistency(qw)
    runExample("ar1_4D", exfolder="../../tmb_examples")
    set.seed(123)
    qw <- checkConsistency(obj, opt$par, n=100)
    qw
}
