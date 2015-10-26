## Copyright (C) 2013-2015 Kasper Kristensen
## License: GPL-2

##' Calculate 1D likelihood profiles wrt. single parameters or more
##' generally, wrt. arbitrary linear combinations of parameters
##' (e.g. contrasts).
##'
##' Given a linear combination \deqn{ t = \sum_{i=1}^n v_i \theta_i } of
##' the parameter vector \eqn{\theta}, this function calculates the
##' likelihood profile of \eqn{t}. By default \eqn{v} is a unit vector
##' determined from \code{name}. Alternatively the linear combination
##' may be given directly (\code{lincomb}).
##'
##' @title Adaptive likelihood profiling.
##' @param obj Object from \code{MakeADFun} that has been optimized.
##' @param name Name or index of a parameter to profile.
##' @param lincomb Optional linear combination of parameters to
##' profile. By default a unit vector corresponding to \code{name}.
##' @param h Initial adaptive stepsize on parameter axis.
##' @param ytol Adjusts the range of the likelihood values.
##' @param ystep Adjusts the resolution of the likelihood profile.
##' @param maxit Max number of iterations for adaptive algorithm.
##' @param slice Do slicing rather than profiling?
##' @param parm.range Valid parameter range.
##' @param trace Trace progress?
##' @param ... Unused
##' @return data.frame with parameter and function values.
##' @seealso \code{\link{plot.tmbprofile}}, \code{\link{confint.tmbprofile}}
##' @examples
##' runExample("simple",thisR=TRUE)
##' ## Parameter names for this model:
##' ## beta   beta   logsdu   logsd0
##'
##' ## Profile wrt. sigma0:
##' prof <- tmbprofile(obj,"logsd0")
##' plot(prof)
##' confint(prof)
##'
##' \dontrun{
##' ## Profile the difference between the beta parameters (name is optional):
##' prof2 <- tmbprofile(obj,name="beta1 - beta2",lincomb = c(1,-1,0,0))
##' plot(prof2)
##' confint(prof2)
##' }
tmbprofile <- function(obj,
                       name,
                       lincomb,
                       h=1e-4,
                       ytol=2,
                       ystep=.1,
                       maxit=ceiling(5*ytol/ystep),
                       parm.range = c(-Inf, Inf),
                       slice=FALSE,
                       trace=TRUE,...){
    ## Cleanup 'obj' when we exit from this function:
    restore.on.exit <- c("last.par.best",
                         "random.start",
                         "value.best",
                         "last.par",
                         "inner.control",
                         "tracemgc")
    oldvars <- sapply(restore.on.exit, get, envir=obj$env, simplify=FALSE)
    restore.oldvars <- function(){
        for(var in names(oldvars)) assign(var, oldvars[[var]], envir=obj$env)
    }
    on.exit(restore.oldvars())
    
    ## Parameter estimate (thetahat)
    par <- obj$env$last.par.best
    if(!is.null(obj$env$random)) par <- par[-obj$env$random]
    
    ## Determine lincomb vector ('lincomb')
    if(missing(lincomb)){
        if (missing(name)) stop("No 'name' or 'lincomb' specified")
        stopifnot(length(name) == 1)
        if(is.numeric(name)){
            lincomb <- as.numeric(1:length(par)==name)
            name <- names(par)[name]
        }
        else if(is.character(name)){
            if (sum(names(par)==name) != 1) stop("'name' is not unique")
            lincomb <- as.numeric(names(par)==name)
        }
        else stop("Invalid name argument")
    } else {
        if (missing(name)) name <- "parameter"
    }
    stopifnot(length(lincomb) == length(par))

    ## Re-parameterize to direction plus (n-1)-dim-subspace
    ##   theta = t*direction + C %*% s
    X <- diag(length(lincomb))
    i <- which(lincomb != 0)[1]
    X[i,] <- lincomb           ## Linear indep. columns
    invX <- solve(X)
    direction <- invX[,i]
    C <- invX[,-i,drop=FALSE]  ## Now t(lincomb) %*% C = 0 !
    that <- sum( lincomb * par )

    ## Start out with initial increment h and ytol.
    ## * Evaluate and store next function value x1=x0+h, y1=f(x1).
    ## * Repeat as long as abs(y1-y.init)<ytol
    ## * If change is too small double the step size h.
    if(slice){ ## Simple slice case
      f <- function(x){
        par <- par + x*direction
        obj$fn(par)
      }
    } else { ## Tough profile case
      f <- function(x){
        par <- par + x*direction
        newfn <- function(par0){
          par <- par + as.vector( C %*% par0 )
          obj$fn(par)
        }
        newgr <- function(par0){
          par <- par + as.vector( C %*% par0 )
          as.vector( obj$gr(par) %*% C )
        }
        ## For inner problem: Use initial guess from previous evaluation
        obj$env$value.best <- Inf
        obj$env$inner.control$trace <- FALSE
        obj$env$tracemgc <- FALSE
        control <- list(step.min=1e-3)
        ans <- nlminb(start,newfn,newgr,control=control)
        start <<- ans$par
        if (trace) cat("Profile value:",ans$objective,"\n")
        ans$objective
      }
    }
    ## Robustify f against failure
    f.original <- f
    f <- function(x){
        y <- try(f.original(x), silent=TRUE)
        if(is(y, "try-error")) y <- NA
        y
    }
    start <- NULL
    evalAlongLine <- function(h){
        start <<- rep(0, length(par)-1)
        x <- 0; y <- f(x)
        if(slice)obj$env$random.start <- expression(last.par[random])
        for(it in 1:maxit){
            yinit <- y[1]
            xcurrent <- tail(x,1)
            ycurrent <- tail(y,1)
            xnext <- xcurrent+h
            if(xnext + that < parm.range[1])                break;
            if(               parm.range[2] < xnext + that) break;
            ynext <- f(xnext)
            x <- c(x,xnext)
            y <- c(y,ynext)
            if( is.na(ynext) )            break;
            if( abs(ynext-yinit) > ytol ) break;
            speedMax <- ystep
            speedMin <-
                if(ynext >= yinit) ystep/4     ## 'tail-part'
                else               ystep/8     ## 'center-part' => slow down
            if( abs(ynext-ycurrent) > speedMax )
                h <- h / 2
            if( abs(ynext-ycurrent) < speedMin )
                h <- h * 2
        }
        ans <- data.frame(x=x+that, y=y)
        names(ans) <- c(name,"value")
        ans
    }
    ans1 <- evalAlongLine(h)
    restore.oldvars()
    ans2 <- evalAlongLine(-h)
    ans <- rbind(ans1,ans2)
    ord <- order(ans[[1]])
    ans <- ans[ord,]
    class(ans) <- c("tmbprofile", class(ans))
    ans
}

##' Plot (negative log) likelihood profile with confidence interval added.
##'
##' @title Plot likelihood profile.
##' @param x Output from \code{\link{tmbprofile}}.
##' @param type Plot type.
##' @param level Add horizontal and vertical lines depicting this confidence level (\code{NULL} disables the lines).
##' @param ... Additional plot arguments.
##' @return NULL
##' @method plot tmbprofile
##' @S3method plot tmbprofile
plot.tmbprofile <- function(x,type="l",level=.95,...){
    plot(as.data.frame(x), type=type, ...)
    if(!is.null(level)){
        hline <- .5*qchisq(level,df=1)
        abline(h=hline+min(x$value), lty="dotted")
        abline(v=confint(x, level=level), lty="dotted")
    }
}

##' Calculate confidence interval from a likelihood profile.
##'
##' @title Profile based confidence intervals.
##' @param object Output from \code{\link{tmbprofile}}.
##' @param parm Not used
##' @param level Confidence level.
##' @param ... Not used
##' @return Lower and upper limit as a matrix.
##' @method confint tmbprofile
##' @S3method confint tmbprofile
confint.tmbprofile <- function (object, parm, level = 0.95, ...){
    i <- which.min(object$value)
    left <- head(object, i)
    right <- tail(object, nrow(object)-i )
    hline <- .5*qchisq(level,df=1) + object$value[i]
    lower <- approx(left[[2]], left[[1]], hline)$y
    upper <- approx(right[[2]], right[[1]], hline)$y
    ans <- t( c(lower=lower, upper=upper) )
    rownames(ans) <- names(object)[1]
    ans
}
