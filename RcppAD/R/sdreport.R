##' After optimization of an AD model, \code{sdreport} is used to calculate standard deviations of
##' all model parameters, including non linear functions of random and fixed effects parameters
##' specified through the ADREPORT() macro from the user template.
##'
##' First, the Hessian wrt. the fixed effect parameter vector (\eqn{\theta}) is calculated.
##' The fixed effects covariance matrix is approximated by
##' \deqn{V(\hat\theta)=\nabla^2 l(\hat\theta)^{-1}}
##' where \eqn{l} denotes the likelihood function (i.e. \code{obj$fn}).
##' 
##' For non-random effect models the standard delta-method is used to calculate the covariance
##' matrix. Let \eqn{\phi(\theta)} denote some non-linear function of \eqn{\theta}. Then
##' \deqn{V(\phi(\hat\theta))\approx \nabla\phi V(\hat\theta) \nabla\phi'}
##' 
##' For random effect models a generalized delta-method is used. First the joint covariance
##' of random and fixed effects is estimated by
##' \deqn{V \pmatrix{ \hat u \cr \hat\theta } \approx
##' \pmatrix{ H_{uu}^{-1} & 0 \cr 0 & 0 } +
##' J V(\hat\theta) J'
##' }
##' where \eqn{H_{uu}} denotes random effect block of the full joint Hessian and \eqn{J}
##' denotes the Jacobian of \eqn{\pmatrix{\hat u(\theta) \cr \theta}} wrt. \eqn{\theta}.
##' Here, the first term represents the expected conditional variance given the fixed effects
##' and the second term represents the variance of the conditional mean wrt. the fixed effects.
##' Now the delta method can be applied on a general non-linear function \eqn{\phi(u,\theta)}
##' of random effects \eqn{u} and fixed effects \eqn{\theta}:
##' \deqn{V(\phi(\hat u,\hat\theta))\approx \nabla\phi V \pmatrix{ \hat u \cr \hat\theta }\nabla\phi'}
##' 
##' @title General sdreport function.
##' @param obj Object returned by \code{MakeADFun}
##' @param par.fixed Optional. Fixed effect parameter estimate (will be known to \code{obj} when an optimization has been carried out).
##' @param hessian.fixed Optional. Hessian wrt. fixed effects (will be calculated from \code{obj} if missing).
##' @return Object of class \code{sdreport}
sdreport <- function(obj,par.fixed=NULL,hessian.fixed=NULL){
  ## Make object to calculate ADREPORT vector
  obj2 <- MakeADFun(obj$env$data,obj$env$parameters,type="ADFun",ADreport=TRUE,DLL=obj$env$DLL)
  r <- obj$env$random
  ## Get full parameter (par), Fixed effects parameter (par.fixed)
  ## and fixed effect gradient (gradient.fixed)
  if(is.null(par.fixed)){ ## Parameter estimate not specified - use best encountered parameter
    par <- obj$env$last.par.best
    if(!is.null(r))par.fixed <- par[-r] else par.fixed <- par
    gradient.fixed <- obj$gr(par.fixed)
  } else {
    gradient.fixed <- obj$gr(par.fixed) ## <-- updates last.par
    par <- obj$env$last.par
  }
  ## Get Hessian wrt. fixed effects (hessian.fixed) and check if positive definite (pdHess).
  if(is.null(hessian.fixed)){
    hessian.fixed <- optimHess(par.fixed,obj$fn,obj$gr) ## Marginal precision of theta.
  }
  pdHess <- !is.character(try(chol(hessian.fixed),silent=TRUE))
  ## ======== Determine case
  ## If no random effects use standard delta method
  simpleCase <- is.null(r)  
  ## Get ADreport vector (phi)
  phi <- try(obj2$fn(par),silent=TRUE)
  if(length(phi)==0){ ## Nothing to report
    simpleCase <- TRUE
  } else { ## Something to report - get derivatives
    Dphi <- obj2$gr(par)
    if(!is.null(r)){
      Dphi.random <- Dphi[,r]
      Dphi.fixed <- Dphi[,-r]
      if(all(Dphi.random==0)){ ## Fall back to simple case
        simpleCase <- TRUE
        Dphi <- Dphi.fixed
      }
    }
  }
  ## ======== Do delta method
  ## Get covariance (cov)
  if(simpleCase){
    cov <- Dphi %*% solve(hessian.fixed) %*% t(Dphi)
  } else {
    hessian.random <- obj$env$spHess(par,random=TRUE)   ## Conditional prec. of u|theta
    L <- obj$env$L.created.by.newton
    if(!is.null(L)){ ## Re-use symbolic factorization if exists
      updateCholesky(L,hessian.random)
      hessian.random@factors <- list(SPdCholesky=L)
    }
    tmp <- solve(hessian.random,t(Dphi.random))
    term1 <- Dphi.random%*%tmp ## first term.
    ## Use columns of tmp as direction for reverse mode sweep
    f <- obj$env$f
    w <- rep(0, length(par))
    reverse.sweep <- function(i){
      w[r] <- tmp[,i]
      -f(par, order = 1, type = "ADGrad",rangeweight = w)[-r]
    }
    A <- t(sapply(seq(length=length(phi)),reverse.sweep)) + Dphi.fixed
    term2 <- A%*%solve(hessian.fixed,t(A)) ## second term
    cov <- term1 + term2
  }
  sd <- sqrt(diag(cov))
  ans <- list(value=phi,sd=sd,cov=cov,par.fixed=par.fixed,
              hessian.fixed=hessian.fixed,pdHess=pdHess,
              gradient.fixed=gradient.fixed)
  class(ans) <- "sdreport"
  ans
}
summary.sdreport <- function(x,...){
  ans1 <- cbind(x$par.fixed,sqrt(diag(solve(x$hessian.fixed))))
  ans2 <- cbind(x$value,x$sd)
  ans <- rbind(ans1,ans2)
  colnames(ans) <- c("Estimate","Std. Error")
  ans
}
print.sdreport <- function(x,...){
  print(summary(x))
  cat("\n")
  if(!x$pdHess){
    cat("Warning:\n")
    cat("Hessian of fixed effects was not positive definite.\n")
  }
  cat("Maximum gradient component:",max(abs(x$gradient.fixed)),"\n")
}
