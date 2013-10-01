##' After optimization of an AD model, \code{sdreport} is used to calculate standard deviations of
##' all model parameters, including non linear functions of random and fixed effects parameters
##' specified through the ADREPORT() macro from the user template.
##'
##' First, the Hessian wrt. the fixed effect parameter vector (\eqn{\theta}) is calculated.
##' The fixed effects covariance matrix is approximated by
##' \deqn{V(\hat\theta)=-\nabla^2 l(\hat\theta)^{-1}}
##' where \eqn{l} denotes the likelihood function (i.e. \code{-obj$fn}).
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
##' where \eqn{H_{uu}} denotes random effect block of the full joint Hessian of \code{obj$fn} and \eqn{J}
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
##' @examples
##' runExample("linreg_parallel",thisR=TRUE) ## Fixed effect example
##' sdreport(obj)
##' runExample("rw",thisR=TRUE)              ## Random effect example
##' rep <- sdreport(obj)
##' summary(rep,"random")                    ## Only random effects
##' summary(rep,"fixed",p.value=TRUE)        ## Only fixed effects
##' summary(rep,"report")                    ## Only report
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
  Vtheta <- solve(hessian.fixed)
  ## Get random effect block of the full joint Hessian (hessian.random) and its
  ## Cholesky factor (L)
  if(!is.null(r)){
    hessian.random <- obj$env$spHess(par,random=TRUE)   ## Conditional prec. of u|theta
    L <- obj$env$L.created.by.newton
    if(!is.null(L)){ ## Re-use symbolic factorization if exists
      updateCholesky(L,hessian.random)
      hessian.random@factors <- list(SPdCholesky=L)
    }
  }
  ## ======== Determine case
  ## If no random effects use standard delta method
  simpleCase <- is.null(r)  
  ## Get ADreport vector (phi)
  phi <- try(obj2$fn(par),silent=TRUE)
  if(is.character(phi) | length(phi)==0){ ## Nothing to report
    simpleCase <- TRUE
    phi <- numeric(0)
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
    if(length(phi)>0){
      cov <- Dphi %*% Vtheta %*% t(Dphi)
    } else cov <- matrix(,0,0)
  } else {
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
    term2 <- A%*%(Vtheta%*%t(A)) ## second term
    cov <- term1 + term2
  }
  ## Output
  sd <- sqrt(diag(cov))
  ans <- list(value=phi,sd=sd,cov=cov,par.fixed=par.fixed,
              cov.fixed=Vtheta,pdHess=pdHess,
              gradient.fixed=gradient.fixed)
  ## ======== Find marginal variances of all random effects i.e. phi(u,theta)=u
  if(!is.null(r)){
    if(is(L,"dCHMsuper")){ ## Required by inverse subset algorithm
      ihessian.random <- .Call("lgc_invQ", L, PACKAGE = "RcppAD")
      diag.term1 <- diag(ihessian.random)
      f <- obj$env$f
      w <- rep(0, length(par))
      reverse.sweep <- function(i){
        w[i] <- 1
        f(par, order = 1, type = "ADGrad",rangeweight = w)[r]
      }
      nonr <- setdiff(seq(length=length(par)),r)
      tmp <- sapply(nonr,reverse.sweep)
      A <- solve(hessian.random,tmp)
      diag.term2 <- rowSums((A %*% Vtheta)*A)
      ans$par.random <- par[r]
      ans$diag.cov.random <- diag.term1 + diag.term2
    } else {
      warning("Could not report sd's of full randomeffect vector.")
    }
  }
  class(ans) <- "sdreport"
  ans
}
summary.sdreport <- function(x,select=c("all","fixed","random","report"),p.value=FALSE,...){
  select <- match.arg(select)
  if(select=="all"){
    fixed <- random <- report <- TRUE
    all <- TRUE
  } else {
    fixed <- random <- report <- FALSE
    assign(select,TRUE)
    all <- FALSE
  }
  ans1 <- ans2 <- ans3 <- NULL
  if(fixed)ans1 <- cbind(x$par.fixed,sqrt(diag(x$cov.fixed)))
  if(random)ans2 <- cbind(x$par.random,sqrt(as.numeric(x$diag.cov.random)))
  if(report)ans3 <- cbind(x$value,x$sd)
  ans <- rbind(ans1,ans2,ans3)
  colnames(ans) <- c("Estimate","Std. Error")
  if(p.value){
    ans <- cbind(ans,p.value=pchisq((ans[,"Estimate"]/ans[,"Std. Error"])^2,df=1,lower=FALSE))
  }
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
