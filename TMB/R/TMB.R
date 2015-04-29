##' Internal TMB Functions
##'
##' Internal TMB functions
##'
##' These are not to be called by the user (or in some cases are just
##' waiting for proper documentation to be written :).
##'
##' @name TMB-internal
##' @aliases checkSparseHessian config dynlib flagsDefaults getUserDLL grepRandomParameters info isParallelTemplate newtonDefaults newtonOption parallelBenchmark plot.parallelBenchmark print.backtrace print.sdreport runSymbolicAnalysis setDefaults sparseHessianFun summary.sdreport tmbOption updateCholesky
##' @rdname TMB-internal
NULL

## Utilities
grepRandomParameters <- function(parameters,random){
  r <- sort(unique(unlist(lapply(random,function(regexp)grep(regexp,names(parameters))))))
  tmp <- lapply(parameters,function(x)x*0)
  tmp[r] <- lapply(tmp[r],function(x)x*0+1)
  which(as.logical(unlist(tmp)))
}

## Guess name of user's loaded DLL code
getUserDLL <- function(){
    dlls <- getLoadedDLLs()
    isTMBdll <- function(dll)!is(try(getNativeSymbolInfo("MakeADFunObject",dll),TRUE),"try-error")
    TMBdll <- sapply(dlls, isTMBdll)
    if(sum(TMBdll) == 0) stop("There are no TMB models loaded (use 'dyn.load').")
    if(sum(TMBdll) >1 ) stop("Multiple TMB models loaded. Failed to guess DLL name.")
    names(dlls[TMBdll])
}

## Update cholesky factorization ( of H+t*I ) avoiding copy overhead
## by writing directly to L(!).
updateCholesky <- function(L,H,t=0){
  .Call("destructive_CHM_update",L,H,as.double(t),PACKAGE="Matrix")
}

##' Construct objective functions with derivatives based on the users c++ template.
##'
##' A call to \code{MakeADFun} will return an object that, based on the users DLL code (specified through \code{DLL}), contains functions to calculate the objective function
##' and its gradient. The object contains the following components:
##' \itemize{
##'   \item \code{par} A default parameter.
##'   \item \code{fn} The likelihood function.
##'   \item \code{gr} The gradient function.
##'   \item \code{report} A function to report all variables reported with the REPORT() macro in the user template.
##'   \item \code{env} Environment with access to all parts of the structure.
##' }
##' and is thus ready for a call to R's \code{optim} function.
##' Data (\code{data}) and parameters (\code{parameters}) are directly read by the user template via the macros beginning with DATA_
##' and PARAMETER_. The order of the PARAMETER_ macros defines the order of parameters in the final objective function.
##' There are no restrictions on the order of random parameters, fixed parameters or data in the template.
##' 
##' Optionally, a simple mechanism for collecting and fixing parameters from R is available through the \code{map} argument. A map is a named list
##' of factors with the following properties:
##' \itemize{
##'   \item names(map) is a subset of names(parameters).
##'   \item For a parameter "p" length(map$p) equals length(parameters$p).
##'   \item Parameter entries with NAs in the factor are fixed.
##'   \item Parameter entries with equal factor level are collected to a common value.
##' }
##' More advanced parameter mapping, such as collecting parameters between different vectors etc., must be implemented from the template.
##' 
##' Random effects are specified via the argument \code{random}: A component of the parameter list is marked as random if its name is matched
##' by any of the characters of the vector \code{random} (Regular expression match is performed if \code{regexp=TRUE}).
##' If some parameters are specified as random effects, these will
##' be integrated out of the objective function via the Laplace approximation. In this situation the functions \code{fn} and \code{gr}
##' automatically perform an optimization of random effects for each function evaluation. This is referred to as
##' the 'inner optimization'. Strategies for choosing initial values of the inner optimization can be controlled
##' via the argument \code{random.start}. The default is \code{expression(last.par.best[random])}
##' where \code{last.par.best} is an internal full parameter vector corresponding to the currently best
##' likelihood. An alternative choice could be \code{expression(last.par[random])} i.e. the random effect optimum of
##' the most recent - not necessarily best - likelihood evaluation. Further control of the inner optimization can
##' be obtained by the argument \code{inner.control} which is a list of control parameters for the inner optimizer
##' \code{newton}. Depending of the inner optimization problem type the following settings are recommended:
##' \enumerate{
##'   \item Quasi-convex: \code{smartsearch=TRUE} (the default).
##'   \item Strictly-convex: \code{smartsearch=FALSE} and \code{maxit=20}.
##'   \item Quadratic: \code{smartsearch=FALSE} and \code{maxit=1}.
##' }
##' 
##' Technically, the user template is processed several times by inserting
##' different types as template parameter, selected by argument \code{type}:
##' \itemize{
##'   \item \code{"ADFun"} Run through the template with AD-types and produce a stack of operations representing the objective function.
##'   \item \code{"Fun"} Run through the template with ordinary double-types.
##'   \item \code{"ADGrad"} Run through the template with nested AD-types and produce a stack of operations representing the objective function gradient.
##' }
##' Each of these are represented by external pointers to c++ structures available in the environment \code{env}.
##'
##' Further objects in the environment \code{env}:
##' \itemize{
##'   \item \code{validpar} Function defining the valid parameter region (by default no restrictions). If an invalid
##' parameter is inserted \code{fn} immediately return NaN.
##'   \item \code{parList} Function to get the full parameter vector of random and fixed effects in a convenient
##' list format.
##'   \item \code{random} An index vector of random effect positions in the full parameter vector.
##'   \item \code{last.par} Full parameter of the latest likelihood evaluation.
##'   \item \code{last.par.best} Full parameter of the best likelihood evaluation.
##'   \item \code{tracepar} Trace every likelihood evaluation ?
##'   \item \code{tracemgc} Trace mgc of every gradient evaluation ?
##'   \item \code{silent} Pass 'silent=TRUE' to all try-calls ?
##' }
##'
##' A high level of tracing information will be output by default when evaluating the objective function and gradient.
##' This is useful while developing a model, but may eventually become annoying. Disable all tracing by passing
##' \code{silent=TRUE} to the \code{MakeADFun} call.
##' 
##' @title Construct objective functions with derivatives based on a compiled c++ template.
##' @param data List of data objects (vectors,matrices,arrays,factors,sparse matrices) required by the user template (Order does not matter and un-used components are allowed).
##' @param parameters List of all parameter objects required by the user template (both random and fixed effects).
##' @param map List defining how to optionally collect and fix parameters - see details.
##' @param type Character vector defining which operation stacks are generated from the users template - see details.
##' @param random Character vector defining the random effect parameters. See also \code{regexp}.
##' @param profile Parameters to profile out of the likelihood (this subset will be appended to \code{random} with Laplace
##' approximation disabled).
##' @param random.start Expression defining the strategy for choosing random effect initial values as function of previous function evaluations - see details.
##' @param hessian Calculate Hessian at optimum?
##' @param method Outer optimization method.
##' @param inner.method Inner optimization method (see function "newton")
##' @param inner.control List controlling inner optimization
##' @param MCcontrol List controlling importance sampler (turned off by default) 
##' @param ADreport Calculate derivatives of macro ADREPORT(vector) instead of objective_function return value?
##' @param atomic Allow tape to contain atomic functions?
##' @param LaplaceNonZeroGradient Allow Taylor expansion around non-stationary point?
##' @param DLL Name of shared object file compiled by user.
##' @param checkParameterOrder Optional check for correct parameter order.
##' @param regexp Match random effects by regular expressions?
##' @param silent Disable all tracing information?
##' @param ... Currently unused.
##' @return List with components (fn,gr, etc) suitable for an optim call.
MakeADFun <- function(data,parameters,map=list(),
                      type=c("ADFun","Fun","ADGrad"[!is.null(random)]),
                      random=NULL,
                      profile=NULL,
                      random.start=expression(last.par.best[random]),
                      hessian=FALSE,method="BFGS",
                      inner.method="newton",
                      inner.control=list(maxit=1000),
                      MCcontrol=list(doMC=FALSE,seed=123,n=100),
                      ADreport=FALSE,
                      atomic=TRUE,
                      LaplaceNonZeroGradient=FALSE, ## Experimental feature: Allow expansion around non-stationary point
                      DLL=getUserDLL(),
                      checkParameterOrder=TRUE, ## Optional check
                      regexp=FALSE,
                      silent=FALSE,
                      ...){
  env <- environment() ## This environment
  if(!is.list(data))
    stop("data must be a list")
  ok <- function(x)(is.matrix(x)|is.vector(x)|is.array(x))&is.numeric(x)
  ok.data <- function(x)ok(x)|is.factor(x)|is(x,"sparseMatrix")|is.list(x)
  check.passed <- function(x){
    y <- attr(x,"check.passed")
    if(is.null(y)) FALSE else y
  }
  if(!check.passed(data)){
    if(!all(sapply(data,ok.data))){
      cat("Problem with these data entries:\n")
      print(which(!sapply(data,ok.data)))
      stop("Only numeric matrices, vectors, arrays, ",
           "factors or lists ",
           "can be interfaced")
    }
  }
  if(!check.passed(parameters)){
    if(!all(sapply(parameters,ok))){
      cat("Problem with these parameter entries:\n")
      print(which(!sapply(parameters,ok)))
      stop("Only numeric matrices, vectors and arrays ",
           "can be interfaced")
    }
  }
  if(length(data)){
    dataSanitize <- function(x){
      if(is.list(x)) return( lapply(x, dataSanitize) )
      if(is(x,"sparseMatrix")){
        x <- as(x,"dgTMatrix")
      } else {
        if(is.factor(x))x <- unclass(x)-1L ## Factors are passed as 0-based integers !!!
        storage.mode(x) <- "double"
      }
      x
    }
    if(!check.passed(data)){
      data <- lapply(data,dataSanitize)
    }
    attr(data,"check.passed") <- TRUE
  }
  if(length(parameters)){
    parameterSanitize <- function(x){
      storage.mode(x) <- "double"
      x
    }
    if(!check.passed(parameters)){
      parameters <- lapply(parameters,parameterSanitize)
    }
    attr(parameters,"check.passed") <- TRUE
  }

  if(checkParameterOrder){
    ## For safety, check that parameter order match the parameter order in user template.
    ## If not, permute parameter list with a warning.
    ## Order in which parameters were requested:
    parNameOrder <- .Call("getParameterOrder",data,parameters,new.env(),PACKAGE=DLL)
    if(!identical(names(parameters),parNameOrder)){
      cat("Order of parameters:\n")
      print(names(parameters))
      cat("Not matching template order:\n")
      print(parNameOrder)
      parameters <- parameters[parNameOrder]
      cat("Your parameter list has been re-ordered.\n(Disable this warning with checkParameterOrder=FALSE)\n")
    }
  }
  
  ## Prepare parameter mapping.
  ## * A parameter map is a factor telling which parameters should be grouped
  ## * NA values are untouched: So user can e.g. set them to zero
  ## * NOTE: CURRENTLY ONLY WORKS ON PARAMETER_ARRAY() !!!
  if(length(map)>0){
    ok <- all(names(map)%in%names(parameters))
    if(!ok)stop("Names in map must correspond to parameter names")
    ok <- all(sapply(map,is.factor))
    if(!ok)stop("map must contain factors")
    ok <- sapply(parameters[names(map)],length)==sapply(map,length)
    if(!all(ok))stop("A map factor length must equal parameter length")
    param.map <- lapply(names(map),
                        function(nam)
                        {
                          ## Shortened parameter
                          ans <- tapply(parameters[[nam]],map[[nam]],mean)
                          if(length(ans)==0)ans <- as.numeric(ans) ## (zero-length case)
                          ## Integer code used to fill short into original shape
                          fnew <- unclass(map[[nam]])
                          fnew[!is.finite(fnew)] <- 0L
                          fnew <- fnew-1L
                          ## Output
                          attr(ans,"shape") <- parameters[[nam]]
                          attr(ans,"map") <- fnew
                          attr(ans,"nlevels") <- length(ans)
                          ans
                        })
    ## Now do the change:
    parameters[names(map)] <- param.map
  }

  ## Utility to get back parameter list in original shape
  parList <- function(x=par[-random],par=last.par){
    ans <- parameters
    nonemp <- sapply(ans,function(x)length(x)>0) ## Workaround utils::relist bug for empty list items
    nonempindex <- which(nonemp)
    skeleton <- as.relistable(ans[nonemp])
    if(any(random)){
      par[-random] <- x
    } else {
      par[] <- x
    }
    li <- relist(par,skeleton)
    reshape <- function(x){
      if(is.null(attr(x,"map")))return(x)
      y <- attr(x,"shape")
      f <- attr(x,"map")
      i <- which(f>=0)
      y[i] <- x[f[i]+1]
      y
    }
    for(i in seq(skeleton)){
      ans[[nonempindex[i]]][] <- as.vector(li[[i]])
    }
    for(i in seq(ans)){
      ans[[i]] <- reshape(ans[[i]])
    }
    ans
  }
  
  #type <- match.arg(type)
  #if("ADFun"%in%type)ptrADFun <- .Call("MakeADFunObject",data,parameters) else ptrADFun <- NULL
  reportenv <- new.env()
  par <- NULL
  last.par.ok <- last.par <- last.par1 <- last.par2 <- last.par.best <- NULL
  value.best <- Inf
  ADFun <- NULL
  Fun <- NULL
  ADGrad <- NULL
  tracepar <- FALSE
  validpar <- function(x)TRUE
  tracemgc <- TRUE
  ## Disable all tracing information
  beSilent <- function(){
      tracemgc <<- FALSE
      inner.control$trace <<- FALSE
      silent <<- TRUE
      cf <- config(DLL=DLL)
      i <- grep("^trace.",names(cf))
      cf[i] <- 0
      cf$DLL <- DLL
      do.call(config, cf)
      NULL
  }
  if(silent)beSilent()

  ## All external pointers are created in function "retape" and can be re-created
  ## by running retape() if e.g. the number of openmp threads is changed.
  retape <- function(){
    if(atomic){ ## FIXME: Then no reason to create ptrFun again later ?
      ## User template contains atomic functions ==>
      ## Have to call "double-template" to trigger tape generation
      Fun <<- .Call("MakeDoubleFunObject",data,parameters,reportenv,PACKAGE=DLL)
      ## Hack: unlist(parameters) only guarantied to be a permutation of the parameter vecter.
      .Call("EvalDoubleFunObject",Fun$ptr,unlist(parameters),control=list(order=as.integer(0)),PACKAGE=DLL)
    }
    if(is.character(profile)){
        random <<- c(random, profile)
    }
    if(is.character(random)){
      if(!regexp){ ## Default: do exact match
        if(!all(random %in% names(parameters))){
          cat("Some 'random' effect names does not match 'parameter' list:\n")
          print(setdiff(random,names(parameters)))
          cat("(Note that regular expression match is disabled by default)\n")
          stop()
        }
        if(any(duplicated(random))){
          cat("Duplicates in 'random' - will be removed\n")
          random <<- unique(random)
        }
        tmp <- lapply(parameters,function(x)x*0)
        tmp[random] <- lapply(tmp[random],function(x)x*0+1)
        random <<- which(as.logical(unlist(tmp)))
      }
      if(regexp){ ## Original regular expression match
        random <<- grepRandomParameters(parameters,random)
        if(length(random)==0){
          cat("Selected random effects did not match any model parameters.\n")
          random <<- NULL
        }
      }
      if(is.character(profile)){
          ## Convert 'profile' to a pointer into random (represented
          ## as logical index vector):
          tmp <- lapply(parameters,function(x)x*0)
          tmp[profile] <- lapply(tmp[profile],function(x)x*0+1)
          profile <<- match( which(as.logical(unlist(tmp))) , random )
          if(any(duplicated(profile))) stop("Profile parameter vector not unique.")
          tmp <- rep(0L, length(random))
          tmp[profile] <- 1L
          profile <<- tmp
      }
      par <<- unlist(parameters)
    }
    if("ADFun"%in%type){
      ADFun <<- .Call("MakeADFunObject",data,parameters,reportenv,
                     control=list(report=as.integer(ADreport)),PACKAGE=DLL)
      par <<- attr(ADFun$ptr,"par")
      last.par <<- par
      last.par1 <<- par
      last.par2 <<- par
      last.par.best <<- par
    }
    if("Fun"%in%type)
      Fun <<- .Call("MakeDoubleFunObject",data,parameters,reportenv,PACKAGE=DLL)
    if("ADGrad"%in%type)
      ADGrad <<- .Call("MakeADGradObject",data,parameters,reportenv,PACKAGE=DLL)
    ## Skip fixed effects from the full hessian ?
    ## * Probably more efficient - especially in terms of memory.
    ## * Only possible if a taped gradient is available - see function "ff" below.
    env$skipFixedEffects <- !is.null(ADGrad)
    delayedAssign("spHess",sparseHessianFun(env, skipFixedEffects=skipFixedEffects ), assign.env = env )
  }

  retape()

  ## Has atomic functions been generated for the tapes ?
  usingAtomics <- function().Call("usingAtomics", PACKAGE=DLL)
  
  f <- function(theta=par,order=0,type=c("ADdouble","double","ADGrad"),
                cols=NULL,rows=NULL,
                sparsitypattern=0,rangecomponent=1,rangeweight=NULL,
                dumpstack=0){
    type <- match.arg(type)
    if(type=="ADdouble"){
      res <- .Call("EvalADFunObject",ADFun$ptr,theta,
                   control=list(
                     order=as.integer(order),
                     hessiancols=as.integer(cols),
                     hessianrows=as.integer(rows),
                     sparsitypattern=as.integer(sparsitypattern),
                     rangecomponent=as.integer(rangecomponent),
                     rangeweight=rangeweight,
                     dumpstack=as.integer(dumpstack)
                     ),
                   PACKAGE=DLL
                   )
      last.par <<- theta
      if(order==1)last.par1 <<- theta
      if(order==2)last.par2 <<- theta
    } else
    if(type=="double"){
      res <- .Call("EvalDoubleFunObject",Fun$ptr,theta,
                   control=list(order=as.integer(order)),PACKAGE=DLL)
    }
    if(type=="ADGrad"){
      res <- .Call("EvalADFunObject",ADGrad$ptr,
                   theta,control=list(order=as.integer(order),
                           hessiancols=as.integer(cols),
                           hessianrows=as.integer(rows),
                           sparsitypattern=as.integer(sparsitypattern),
                           rangecomponent=as.integer(rangecomponent),
                           rangeweight=rangeweight,
                           dumpstack=as.integer(dumpstack)),PACKAGE=DLL)
    }
    res
  }

  h <- function(theta=par,order=0,hessian,L,...){
    if(order==0){
      ##logdetH <- determinant(hessian)$mod
      logdetH <- 2*determinant(L)$mod
      ans <- f(theta,order=0)+
        .5*logdetH - length(random)/2*log(2*pi)
      if(LaplaceNonZeroGradient){
        grad <- f(theta,order=1)[random]
        ans <- ans - .5*sum(grad*as.numeric(solve(L,grad)))
      }
    }
    if(order==1){
      if(LaplaceNonZeroGradient)stop("Not correct for LaplaceNonZeroGradient=TRUE")
      ##browser()
      e <- environment(spHess)
      solveSubset <- function(L).Call("tmb_invQ",L,PACKAGE="TMB")
      solveSubset2 <- function(L).Call("tmb_invQ_tril_halfdiag",L,PACKAGE="TMB")
      ## FIXME: The following two lines are not efficient:
      ## 1. ihessian <- tril(solveSubset(L))
      ## 2. diag(ihessian) <- .5*diag(ihessian)
      ## Make option to solveSubset to return lower triangular part
      ## with diagonal halved. As it is now the output of solveSubset is
      ## symm _with upper storage_ (!) (side effect of cholmod_ptranspose)
      ## therefore tril takes long time. Further, "diag<-" is too slow.
      ## FIXED! :
      ihessian <- solveSubset2(L)
      ## Profile case correction (1st order case only)
      if(!is.null(profile)){
          ## Naive way:
          ##   ihessian[profile,] <- 0
          ##   ihessian[,profile] <- 0
          ## However, this would modify sparseness pattern and also not
          ## account for 'ihessian' being permuted:
          perm <- L@perm+1L
          ihessian <- .Call("tmb_sparse_izamd", ihessian, profile[perm], 0.0, PACKAGE="TMB")
      }
      
      ## General function to lookup entries A subset B.
      ## lookup.old <- function(A,B){
      ##   A <- as(tril(A),"dtTMatrix")
      ##   B <- as(tril(B),"dtTMatrix")
      ##   match(paste(A@i,A@j),paste(B@i,B@j))
      ## }
      ## General function to lookup entries A in B[r,r] assuming pattern of A
      ## is subset of pattern of B[r,r].
      lookup <- function(A,B,r=NULL){
        A <- tril(A);B <- tril(B)
        B@x[] <- seq.int(length.out=length(B@x)) ## Pointers to full B matrix (FIXME: what if length(B@x)>2^32 ? )
        B <- forceSymmetric(B)
        if(!is.null(r))B <- B[r,r,drop=FALSE] ## Reduce to have same dim as A
        m <- .Call("match_pattern",A,B,PACKAGE="TMB") ## Same length as A@x with pointers to B@x
        B@x[m]
      }
      if(is.null(e$ind1)){
        ## hessian: Hessian of random effect part only.
        ## ihessian: Inverse subset of hessian (same dim but larger pattern!).
        ## Hfull: Pattern of full hessian including fixed effects.
        if (!silent) cat("Matching hessian patterns... ")
        iperm <- Matrix::invPerm(L@perm+1L)
        e$ind1 <- lookup(hessian,ihessian,iperm) ## Same dimensions
        e$ind2 <- lookup(hessian,e$Hfull,random)  ## Note: dim(Hfull)>dim(hessian) !
        if (!silent) cat("Done\n")
      }
      w <- rep(0,length=length(e$Hfull@x))
      w[e$ind2] <- ihessian@x[e$ind1]
      ## Reverse mode evaluate ptr in rangedirection w
      ## now gives .5*tr(Hdot*Hinv) !!
      ans <- as.vector(f(theta,order=1))+
        .Call("EvalADFunObject",e$ADHess$ptr,theta,
                   control=list(
                     order=as.integer(1),
                     hessiancols=as.integer(0),
                     hessianrows=as.integer(0),
                     sparsitypattern=as.integer(0),
                     rangecomponent=as.integer(1),
                     rangeweight=as.double(w),
                     dumpstack=as.integer(0)
                     ),PACKAGE=DLL
                   )

      
    }
    return(ans)
  }

  ff <- function(par.fixed=par[-random],order=0,...){
    names(par.fixed) <- names(par[-random])
    f0 <- function(par.random,order=0,...){
      par[random] <- par.random
      par[-random] <- par.fixed
      res <- f(par,order=order,...)
      switch(order+1,res,res[random],res[random,random])
    }
    ## sparse hessian
    H0 <- function(par.random){
      par[random] <- par.random
      par[-random] <- par.fixed
      #spHess(par)[random,random,drop=FALSE]
      spHess(par,random=TRUE)
    }
    if(inner.method=="newton"){
      #opt <- newton(eval(random.start),fn=f0,gr=function(x)f0(x,order=1),
      #              he=function(x)f0(x,order=2))
      opt <- try( do.call("newton",c(list(par=eval(random.start),
                                      fn=f0,
                                      gr=function(x)f0(x,order=1),
                                      ##he=function(x)f0(x,order=2)),
                                      he=H0,env=env),
                                 inner.control)
                          ), silent=silent
                 )
      if(is.character(opt))return(NaN)
    } else{  
      opt <- optim(eval(random.start),fn=f0,gr=function(x)f0(x,order=1),
                   method=inner.method,control=inner.control)
    }
    par[random] <- opt$par
    par[-random] <- par.fixed

    ## HERE! - update hessian and cholesky
    if(!skipFixedEffects){ ## old way
      hess <- spHess(par) ## Full hessian
      hessian <- hess[random,random] ## Subset
    } else {
      hessian <- spHess(par,random=TRUE)
    }
    ## Profile case correction (0 and 1st order)
    if( !is.null(profile) ){
        ## Naive way:
        ##   hessian[profile, ] <- 0
        ##   hessian[, profile] <- 0
        ##   diag(hessian)[profile] <- 1
        ## However, this would modify sparseness pattern:
        hessian <- .Call("tmb_sparse_izamd", hessian, profile, 1.0, PACKAGE="TMB")
    }
    ## Update Cholesky:
    if(inherits(env$L.created.by.newton,"dCHMsuper")){
      L <- env$L.created.by.newton
      ##.Call("destructive_CHM_update",L,hessian,as.double(0),PACKAGE="Matrix")
      updateCholesky(L,hessian)
    } else
    L <- Cholesky(hessian,perm=TRUE,LDL=FALSE,super=TRUE)

    
    if(order==0){
      res <- h(par,order=0,hessian=hessian,L=L)
      ## Profile case correction
      if(!is.null(profile)){
          res <- res + sum(profile)/2*log(2*pi)
      }
      if(is.finite(res)){
        if(res<value.best){
          last.par.best <<- par; value.best <<- res
        }
      }
    }
    if(order==1){
      #hess <- f(par,order=2,cols=random)
      #hess <- spHess(par)##[,random,drop=FALSE]
      grad <- h(par,order=1,hessian=hessian,L=L)
      #res <- grad[-random] - t(grad[random])%*%solve(hess[random,random])%*%hess[random,-random]
      #res <- grad[-random] - t(grad[random])%*%solve(hess[random,])%*%t(hess[-random,])

      ## Profile case correction. The remaining calculations must be
      ## done with the original hessian (which has been destroyed)
      if(!is.null(profile)){
          ## Update hessian and Cholesky:
          if(!skipFixedEffects){ ## old way
              hess <- spHess(par) ## Full hessian
              hessian <- hess[random,random] ## Subset
          } else {
              hessian <- spHess(par,random=TRUE)
          }
          updateCholesky(L,hessian)
      }
      
      ## res <- grad[-random] -
      ##   hess[-random,random]%*%as.vector(solve(hess[random,random],grad[random]))

      if(!skipFixedEffects){
        ## Relies on "hess[-random,random]" !!!!!
        res <- grad[-random] -
          hess[-random,random]%*%as.vector(solve(L,grad[random]))        
      } else {
        ## Smarter: Do a reverse sweep of ptrADGrad
        w <- rep(0,length(par))
        w[random] <- as.vector(solve(L,grad[random]))
        res <- grad[-random] -
          f(par,order=1,type="ADGrad",rangeweight=w)[-random]
      }
      
      res <- drop(res)
    }
    if(order==2){
      n <- length(par); nr <- length(random); nf <- n-nr
      fixed <- setdiff(1:n,random)
      D1h <- h(par,order=1) ## 1*n
      D2h <- h(par,order=2) ## n*n
      D2f <- f(par,order=2,cols=random) ## n*nr
      D3f <- sapply(random,function(i)
                    f(par,type="ADGrad",
                      order=2,rangecomponent=i))  ## n^2 * nr
      D1eta <- -t(D2f[-random,]%*%solve(D2f[random,]))  ## nr*nf
      D3f.D1eta <- D3f%*%D1eta ## n^2 * nf
      dim(D3f.D1eta) <- c(n,n,nf)
      dim(D3f) <- c(n,n,nr)
      D3f.fixed <- D3f[fixed,,] ##nf*n*nr
      D2eta <- sapply(1:nf,function(i){
        -solve(D2f[random,]) %*%
        ( t(D3f.fixed[i,fixed,]) +  D3f.D1eta[random,fixed,i] +
        ( D3f.fixed[i,random,] + D3f.D1eta[random,random,i] ) %*% D1eta )
        }) # nr*nf*nf
      dim(D2eta) <- c(nr,nf,nf)
      D2h.fixed <- D2h[fixed,] #nf*n
      res <- sapply(1:nf,function(i){
        D2h.fixed[i,fixed] + t( D2h.fixed[,random] %*% D1eta[,i] ) +
          ( t(D2h.fixed[i,random]) +
           t(D2h[random,random] %*% D1eta[,i]) ) %*% D1eta +
             D1h[,random] %*% D2eta[,,i]
      })
      attr(res,"D2eta") <- D2eta
      attr(res,"D1eta") <- D1eta
      #attr(res,"D1h") <- D1h
      #attr(res,"D2h") <- D2h
      #attr(res,"D2f") <- D2f
      #attr(res,"D3f") <- D3f
    }
    if(all(is.finite(res)))last.par.ok <<- par
    return(res)
  }

  ## Monte Carlo improvement of Laplace approximation
  ## Importance sampling from *fixed* reference measure determined
  ## by parameter "par0". Assumptions:
  ## * We know how to sample from the measure - see rmvnorm.
  ## * We know how to evaluate the density of the samples - see logdmvnorm.
  ## * Eventually "par0" will stabilize and become independent of the fixed
  ##   effects "par", so that derivatives of the sample density and the samples
  ##   are zero wrt. the fixed effects.
  MC <- function(par=last.par,       ## Point in which we are evaluating the likelihood
                 par0=last.par.best, ## Parameter for proposal distribution
                 n=100,              ## Number of samples
                 order=0,            ## Derivative order
                 seed=NULL,          ## Random seed
                 antithetic=TRUE,    ## Reduce variance
                 keep=FALSE,         ## Keep samples and fct evals
                 phi=NULL,           ## Function to calculate mean of
                 ...){
    if(is.numeric(seed))set.seed(seed)
    ## Clean up on exit
    last.par.old <- last.par
    last.par.best.old <- last.par.best
    on.exit({last.par <<- last.par.old;
             last.par.best <<- last.par.best.old})
    ## Update Cholesky needed by reference measure
    h <- spHess(par0,random=TRUE)
    L <- L.created.by.newton
    updateCholesky(L,h)             ## P %*% h %*% Pt = L %*% Lt
    rmvnorm <- function(n){
        u <- matrix(rnorm(ncol(L)*n),ncol(L),n)
        u <- solve(L,u,system="Lt") ## Solve Lt^-1 %*% u
        u <- solve(L,u,system="Pt") ## Multiply Pt %*% u
        as.matrix(u)
    }
    logdmvnorm <- function(u){
        logdetH <- 2*determinant(L,logarithm=TRUE)$modulus
        ans <- nrow(h)*log(1/sqrt(2*pi))+.5*logdetH-.5*colSums(u*as.matrix(h%*%u))
        ans
    }
    eval.target <- function(u,order=0){
      par[random] <- u
      f(par,order=order)
    }
    samples <- rmvnorm(n)
    if(antithetic)samples <- cbind(samples,-samples) ## Antithetic variates
    log.density.propose <- logdmvnorm(samples)
    samples <- samples+par0[random]
    log.density.target <- -apply(samples,2,eval.target)
    log.density.target[is.nan(log.density.target)] <- -Inf
    I <- log.density.target - log.density.propose
    M <- max(I)
    if(order>=1){
      vec <- exp(I-M)
      p <- vec/sum(vec)
      i <- (p>0)
      p <- p[i]
      I1 <- apply(samples[,i,drop=FALSE],2,eval.target,order=1)[-random,,drop=FALSE]
      gr <- as.vector(I1 %*% p)
      if(order==1)return(gr)
      ## I1I1 <- t(apply(I1,1,function(x)x%*%t(x)))
      ## I2 <- t(apply(samples,1,function(x)eval.target(x,order=2)[-random,-random]))
      ## h <- colMeans(vec*(-I1I1+I2))/mean(vec)+as.vector(gr)%*%t(as.vector(gr))
      ## if(order==2)return(h)
    }
    if(!is.null(phi)){
      phival <- apply(samples,2,phi)
      if(is.null(dim(phival)))phival <- t(phival)
      p <- exp(I-M); p <- p/sum(p)
      ans <- phival %*% p
      return(ans)
    }
    value <- -log(mean(exp(I-M)))-M
    ci <- 1.96*sd(exp(I-M))/sqrt(n)
    attr(value,"confint") <- -log(mean(exp(I-M))+c(lower=ci,upper=-ci))-M
    if(keep){
        attr(value,"samples") <- samples
        attr(value,"nlratio") <- -I
    }
    value
  }

  report <- function(par=last.par){
    f(par,order=0,type="double")
    as.list(reportenv)
  }

  if(is.null(random)){  ## Output if pure fixed effect model
    return(list(par=par,
                fn=function(x=last.par,...){
                  if(tracepar){cat("par:\n");print(x)}
                  if(!validpar(x))return(NaN)
                  res <- f(x,order=0)
                  if(!ADreport){
                    if(is.finite(res)){
                      if(res<value.best){
                        last.par.best <<- x; value.best <<- res
                      }
                    }
                  }
                  res
                },
                gr=function(x=last.par,...){
                  ans <- f(x,order=1)
                  if(tracemgc)cat("outer mgc: ",max(abs(ans)),"\n")
                  ans
                },
                he=function(x=last.par,atomic=usingAtomics()){
                    ## If no atomics on tape we have all orders implemented:
                    if(!atomic) return( f(x,order=2) )
                    ## Otherwise, get Hessian as 1st order derivative of gradient:
                    if(is.null(ADGrad))
                        ADGrad <<- .Call("MakeADGradObject",data,parameters,reportenv,PACKAGE=DLL)
                    f(x,type="ADGrad",order=1)
                },
                hessian=hessian,method=method,
                retape=retape,env=env,
                report=report,...))
  }
  if(!is.null(random)){  ## Output if random effect model
    return(list(par=par[-random],
                fn=function(x=last.par[-random],...){
                  if(tracepar){cat("par:\n");print(x)}
                  if(!validpar(x))return(NaN)
                  ans <- try({
                    if(MCcontrol$doMC){
                      ff(x,order=0)
                      MC(last.par,n=MCcontrol$n,seed=MCcontrol$seed,order=0)
                    } else
                    ff(x,order=0)
                  },silent=silent)
                  if(is.character(ans))NaN else ans
                },
                gr=function(x=last.par[-random],...){
                  ans <- {
                    if(MCcontrol$doMC){
                      ff(x,order=0)
                      MC(last.par,n=MCcontrol$n,seed=MCcontrol$seed,order=1)
                    } else
                    ff(x,order=1)
                  }
                  if(tracemgc)cat("outer mgc: ",max(abs(ans)),"\n")
                  ans
                },
                he=function(x=last.par[-random],...){
                  stop("Hessian not yet implemented for models with random effects.")
                  if(MCcontrol$doMC){
                    ff(x,order=0)
                    MC(last.par,n=MCcontrol$n,seed=MCcontrol$seed,order=2)
                  } else
                  ff(x,order=2)
                },
                hessian=hessian,method=method,
                retape=retape,env=env,
                report=report,...))
  }
   
  return(env)
}

.removeComments <- function(x){
  x <- paste(x,collapse="\n")
  remlong <- function(x)gsub("/\\*.*?\\*/","",x)
  remshort <- function(x)gsub("//[^\n]*\n","\n",x)
  x <- remshort(remlong(x))
  strsplit(x,"\n")[[1]]
}

isParallelTemplate <- function(file){
  code <- readLines(file)
  code <- .removeComments(code)
  length(grep("^[ \t]*PARALLEL_",code))>0  ||
  length(grep("^[ \t]*parallel_accumulator",code))>0
}

##' Control number of openmp threads.
##'
##' @title Control number of openmp threads.
##' @param n Requested number of threads, or \code{NULL} to just read the current value.
##' @return Number of threads.
openmp <- function(n=NULL){
  if(!is.null(n))n <- as.integer(n)
  .Call("omp_num_threads",n,PACKAGE="TMB")
}

##' Compile a c++ template into a shared object file. OpenMP flag is set if the template is detected to be parallel.
##'
##' TMB relies on R's built in functionality to create shared libraries independent on the platform.
##' A template is compiled by \code{compile("template.cpp")}, which will call R's makefile with appropriate
##' preprocessor flags.
##' Compiler and compiler flags can be stored in a configuration file. In order of precedence either via
##' the file pointed at by R_MAKEVARS_USER or the file ~/.R/Makevars if it exists.
##' Additional configuration variables can be set with \code{...} argument, which will overwrite any
##' previous selections.
##' @title Compile a c++ template to DLL suitable for MakeADFun.
##' @param file c++ file.
##' @param flags Character with compile flags.
##' @param safebounds Turn on preprocessor flag for bound checking?
##' @param safeunload Turn on preprocessor flag for safe DLL unloading?
##' @param openmp Turn on openmp flag? Auto detected for parallel templates.
##' @param libtmb Use precompiled TMB library if available (to speed up compilation)?
##' @param ... Passed as Makeconf variables.
compile <- function(file,flags="",safebounds=TRUE,safeunload=TRUE,
                    openmp=isParallelTemplate(file),libtmb=TRUE,...){
  if(.Platform$OS.type=="windows"){
    ## Overload system.file
    system.file <- function(...){
      ans <- base::system.file(...)
      ans <- chartr("\\", "/", shortPathName(ans))
      ans
    }
  }
  ## libtmb existence
  if(!openmp){
    libtmb <- libtmb && file.exists(system.file(dynlib("libs/libTMB"),package="TMB"))
  }
  if(openmp){
    libtmb <- libtmb && file.exists(system.file(dynlib("libs/libTMBomp"),package="TMB"))
  }
  ## Function to create temporary makevars, Note:
  ## * R_MAKEVARS_USER overrules all other Makevars in tools:::.shlib_internal
  oldmvuser <- mvuser <- Sys.getenv("R_MAKEVARS_USER",NA)
  if(is.na(oldmvuser)){
    on.exit(Sys.unsetenv("R_MAKEVARS_USER"))
  } else {
    on.exit(Sys.setenv(R_MAKEVARS_USER=oldmvuser))
  }
  if(is.na(mvuser) && (file.exists(f <- path.expand("~/.R/Makevars"))))mvuser <- f
  if(!is.na(mvuser)){
    cat("Note: Using Makevars in",mvuser,"\n")
  }
  makevars <- function(...){
    file <- tempfile()
    args <- unlist(list(...))
    txt <- paste(names(args),args,sep="=")
    if(!is.na(mvuser)){
      if(file.exists(mvuser)){
        txt <- c(readLines(mvuser),txt)
      }
    }
    writeLines(txt,file)
    Sys.setenv(R_MAKEVARS_USER=file)
    file
  }
  ## Check that libname is valid C entry.
  libname <- sub("\\.[^\\.]*$","",basename(file))
  if(safeunload){
    valid <- c(letters[1:26],LETTERS[1:26],0:9,"_")
    invalid <- setdiff(unique(strsplit(libname,"")[[1]]),valid)
    if(length(invalid)>0){
      cat("Your library name has invalid characters:\n")
      print(invalid)
      cat("It is recommended to replace invalid characters by underscore.\n")
      cat("Alternatively compile with safeunload=FALSE (not recommended).\n")
      stop()
    }
  }
  ## On windows the DLL must be unloaded before compiling
  if(.Platform$OS.type=="windows"){
    tr <- try(dyn.unload(dynlib(libname)),silent=TRUE)
    if(!is(tr,"try-error"))cat("Note: Library",paste0("'",dynlib(libname),"'"),"was unloaded.\n")
  }
  ## Includes and preprocessor flags specific for the template
  ppflags <- paste(paste0("-I",system.file("include",package="TMB")),
                   "-DTMB_SAFEBOUNDS"[safebounds],
                   paste0("-DLIB_UNLOAD=R_unload_",libname)[safeunload],
                   "-DWITH_LIBTMB"[libtmb]
                   )
  ## Makevars specific for template
  mvfile <- makevars(PKG_CPPFLAGS=ppflags,
                     PKG_LIBS=paste(
                       "$(SHLIB_OPENMP_CXXFLAGS)"[openmp],
                       system.file(dynlib("libs/libTMB"),package="TMB")[libtmb && !openmp],
                       system.file(dynlib("libs/libTMBomp"),package="TMB")[libtmb && openmp] ),
                     PKG_CXXFLAGS="$(SHLIB_OPENMP_CXXFLAGS)"[openmp],
                     CXXFLAGS=flags[flags!=""], ## Optionally overwrite cxxflags
                     ...
                     )
  on.exit(file.remove(mvfile),add=TRUE)
  status <- tools:::.shlib_internal(file)
  if(status!=0)stop("Compilation failed")
  status
}

##' Precompile the TMB library
##'
##' The precompilation should only be run once, typically right after installation of TMB.
##' Note that the precompilation requires write access to the TMB package folder.
##' Two versions of the library - with/without the openmp flag - will be generated. After this,
##' compilation times of templates should be reduced.
##' \itemize{
##' \item To precompile on Linux run \code{precompile()}.
##' \item To precompile on OS X run \code{precompile(PKG_LIBS = "-install_name `pwd`/$@@")}.
##' }
##' Note that precompilation has side effects: It is not possible to work with more than one
##' model at a time for a single R instance.
##' @title Precompile the TMB library in order to speed up compilation of templates.
##' @param ... Passed to \code{compile}.
precompile <- function(...){
  owdir <- getwd()
  on.exit(setwd(owdir))
  folder <- system.file("libs",package="TMB")
  setwd(folder)
  writeLines("#include <TMB.hpp>","libTMB.cpp")
  cat("Compiling serial version\n")
  compile("libTMB.cpp",safeunload=FALSE,libtmb=FALSE,...)
  file.remove("libTMB.cpp")
  writeLines("#include <TMB.hpp>","libTMBomp.cpp")
  cat("Compiling parallel version\n")
  compile("libTMBomp.cpp",openmp=TRUE,safeunload=FALSE,libtmb=FALSE,...)
  file.remove("libTMBomp.cpp")
}

## Add dynlib extension
dynlib <- function(x)paste0(x,.Platform$dynlib.ext)

##' Create a cpp template to get started.
##'
##' This function generates a c++ template with a header and include statement. Here is a brief
##' overview of the c++ syntax used to code the objective function.
##' 
##' Macros to read data and declare parameters:
##'  \tabular{lll}{
##'     \bold{Template Syntax}    \tab     \bold{C++ type}            \tab    \bold{R type} \cr
##'     DATA_VECTOR(name)         \tab     vector<Type>               \tab    vector        \cr
##'     DATA_MATRIX(name)         \tab     matrix<Type>               \tab    matrix        \cr
##'     DATA_SCALAR(name)         \tab     Type                       \tab    numeric(1)    \cr
##'     DATA_INTEGER(name)        \tab     int                        \tab    integer(1)    \cr
##'     DATA_FACTOR(name)         \tab     vector<int>                \tab    factor        \cr
##'     DATA_IVECTOR(name)        \tab     vector<int>                \tab    integer       \cr
##'     DATA_SPARSE_MATRIX(name)  \tab     Eigen::SparseMatrix<Type>  \tab    dgTMatrix     \cr
##'     DATA_ARRAY(name)          \tab     array<Type>                \tab    array         \cr
##'     PARAMETER_MATRIX(name)    \tab     matrix<Type>               \tab    matrix        \cr
##'     PARAMETER_VECTOR(name)    \tab     vector<Type>               \tab    vector        \cr
##'     PARAMETER_ARRAY(name)     \tab     array<Type>                \tab    array         \cr
##'     PARAMETER(name)           \tab     Type                       \tab    numeric(1)    \cr
##'  }
##'
##' Basic calculations:
##'  \tabular{ll}{
##'     \bold{Template Syntax}    \tab   \bold{Explanation}                     \cr
##'     REPORT(x)                 \tab   Report x back to R                     \cr
##'     ADREPORT(x)               \tab   Report x back to R with derivatives    \cr
##'     vector<Type> v(n1);       \tab   R equivalent of v=numeric(n1)          \cr
##'     matrix<Type> m(n1,n2);    \tab   R equivalent of m=matrix(0,n1,n2)      \cr
##'     array<Type> a(n1,n2,n3);  \tab   R equivalent of a=array(0,c(n1,n2,n3)) \cr
##'     v+v,v-v,v*v,v/v           \tab   Pointwise binary operations            \cr
##'     m*v                       \tab   Matrix-vector multiply                 \cr
##'     a.col(i)                  \tab   R equivalent of a[,,i]                 \cr
##'     a.col(i).col(j)           \tab   R equivalent of a[,j,i]                \cr
##'     a(i,j,k)                  \tab   R equivalent of a[i,j,k]               \cr
##'     exp(v)                    \tab   Pointwise math                         \cr
##'     m(i,j)                    \tab   R equivalent of m[i,j]                 \cr
##'     v.sum()                   \tab   R equivalent of sum(v)                 \cr
##'     m.transpose()             \tab   R equivalent of t(m)                   \cr
##'  }
##'
##' Some distributions are avaliable as c++ templates with syntax close to R's distributions:
##' \tabular{ll}{
##'    \bold{Function header}                \tab \bold{Distribution}                      \cr
##'    dnbinom2(x,mu,var,int give_log=0)     \tab Negative binomial with mean and variance \cr
##'    dpois(x,lambda,int give_log=0)        \tab Poisson distribution as in R             \cr
##'    dlgamma(y,shape,scale,int give_log=0) \tab log-gamma distribution                   \cr
##'    dnorm(x,mean,sd,int give_log=0)       \tab Normal distribution as in R              \cr
##' }
##' @title Create cpp template to get started.
##' @param file Optional name of cpp file.
##' @examples
##' template()
template <- function(file=NULL){
  x <- readLines(system.file("template.cpp",package="TMB"))
  if(!is.null(file)){
    if(file.exists(file))stop("File '",file,"' exists")
    writeLines(x,file)
  }
  else cat(paste(x,collapse="\n"))
}

##' Create a skeleton of required R-code once the cpp template is ready.
##'
##' @title Create minimal R-code corresponding to a cpp template.
##' @param file cpp template file.
##' @examples
##' file <- system.file("examples/simple.cpp", package = "TMB")
##' Rinterface(file)
Rinterface <- function(file){
  libname <- sub("\\.[^\\.]*$", "", basename(file))
  x <- readLines(file)
  x <- .removeComments(x)
  items2list <- function(items){
    if(length(items)==0)return("list(),")
    paste0("list(\n",paste(paste0("  ",items,"=  "),collapse=",\n"),"\n ),")
  }
  ## Data
  dataregexp <- "^[ ]*DATA_.*?\\((.*?)\\).*"
  datalines <- grep(dataregexp,x,value=TRUE)
  dataitems <- sub(dataregexp,"\\1",datalines)
  ## Parameters
  parameterregexp <- "^[ ]*PARAMETER.*?\\((.*?)\\).*"
  parameterlines <- grep(parameterregexp,x,value=TRUE)
  parameteritems <- sub(parameterregexp,"\\1",parameterlines)
  libname <- paste0("\"",libname,"\"")
  txt <- c("library(TMB)",
           paste0("dyn.load(dynlib(",libname,"))"),
           "MakeADFun(",
           paste0(" data=",items2list(dataitems)),
           paste0(" parameters=",items2list(parameteritems)),
           paste0(" DLL=",libname),
           ")\n"
           )
  cat(paste(txt,collapse="\n"))
}

## Get som info about the ADFun pointers
info <- function(obj){
  if(!is.environment(obj$env))stop("Wrong object")
  env <- obj$env
  fun <- function(name){
    if(is.null(get(name,env)))return(NA)
    unlist(.Call("InfoADFunObject",get(name,env),PACKAGE=obj$env$DLL))
  }
  names <- ls(env,pattern="ptr")
  lapply(names,fun)
}

## Recommended settings:
## * General non-convex case: smartsearch=TRUE
## * Strictly convex case:    smartsearch=FALSE and maxit=20
## * Quadratic case:          smartsearch=FALSE and maxit=1


##' Generalized newton optimizer used for the inner optimization problem.
##'
##' If \code{smartsearch=FALSE} this function performs an ordinary newton optimization
##' on the function \code{fn} using an exact sparse hessian function.
##' A fixed stepsize may be controlled by \code{alpha} so that the iterations are
##' given by:
##' \deqn{u_{n+1} = u_n - \alpha f''(u_n)^{-1}f'(u_n)}
##'
##' If \code{smartsearch=TRUE} the hessian is allowed to become negative definite
##' preventing ordinary newton iterations. In this situation the newton iterations are performed on
##' a modified objective function defined by adding a quadratic penalty around the expansion point \eqn{u_0}:
##' \deqn{f_{t}(u) = f(u) + \frac{t}{2} \|u-u_0\|^2}{f_t(u) = f(u) + t/2 |u-u_0|^2}
##' This functions hessian ( \eqn{f''(u)+t I} ) is positive definite for \eqn{t} sufficiently
##' large. The value \eqn{t} is updated at every iteration: If the hessian is positive definite \eqn{t} is
##' decreased, otherwise increased. Detailed control of the update process can be obtained with the
##' arguments \code{ustep}, \code{power} and \code{u0}.
##' @title Generalized newton optimizer.
##' @param par Initial parameter.
##' @param fn Objective function.
##' @param gr Gradient function.
##' @param he Sparse hessian function.
##' @param trace Print tracing information?
##' @param maxit Maximum number of iterations.
##' @param tol Convergence tolerance.
##' @param alpha Newton stepsize in the fixed stepsize case.
##' @param smartsearch Turn on adaptive stepsize algorithm for non-convex problems?
##' @param mgcmax Refuse to optimize if the gradient is too steep.
##' @param super Supernodal Cholesky?
##' @param silent Be silent?
##' @param ustep Adaptive stepsize initial guess between 0 and 1.
##' @param power Parameter controlling adaptive stepsize.
##' @param u0 Parameter controlling adaptive stepsize.
##' @param grad.tol Gradient convergence tolerance. 
##' @param step.tol Stepsize convergence tolerance.
##' @param tol10 Try to exit if last 10 iterations not improved more than this.
##' @param env Environment for cached Cholesky factor.
##' @param ... Currently unused.
##' @return List with solution similar to \code{optim} output.
newton <- function (par,fn,gr,he,
                    trace = newtonOption("trace"),
                    maxit = newtonOption("maxit"),
                    tol=newtonOption("tol"),
                    alpha=1,
                    smartsearch=newtonOption("smartsearch"),
                    mgcmax=newtonOption("mgcmax"),
                    super=TRUE,
                    silent=TRUE,
                    ustep = 1, ## Start out optimistic: Newton step
                    power=.5, ## decrease=function(u)const*u^power
                    u0=1e-4,  ## Increase u=0 to this value  
                    grad.tol=tol,
                    step.tol=tol,
                    tol10=1e-3, ## Try to exit if last 10 iterations not improved much
                    env=environment(),
                    ...)
{
  ## Test if a Cholesky factor is present inside the environment of "he" function.
  ## If not - create one...
  if(is.null(env$L.created.by.newton)){
    h.pattern <- he(par)
    ## Make sure Cholesky is succesful
    h.pattern@x[] <- 0
    diag(h.pattern) <- 1
    env$L.created.by.newton <- Cholesky(h.pattern,super=super)
  }
  L <- env$L.created.by.newton
  chol.solve <- function(h,g){
    ##.Call("destructive_CHM_update",L,h,as.double(0),PACKAGE="Matrix")
    updateCholesky(L,h)
    as.vector(solve(L,g))
  }
  optimize <- stats::optimize
  nam <- names(par)
  par <- as.vector(par)
  g <- h <- NULL
  ## pd.check: Quick test for hessian being positive definite
  iterate <- function(par,pd.check=FALSE) {
    if(file.exists(".Rbreakpoint"))browser() ## secret backdoor to poke around...
    if(pd.check){
      if(is.null(h))return(TRUE)
      h <<- he(par) ## Make sure hessian is updated
      tmp <- try( updateCholesky(L,h) , silent=silent)
      return( !inherits(tmp,"try-error") )
    }
    g <<- as.vector(gr(par))
    if(any( !is.finite(g) ))stop("Newton dropout because inner gradient had non-finite components.")
    if(is.finite(mgcmax))
      if(max(abs(g))>mgcmax)stop("Newton dropout because inner gradient too steep.")
    if(max(abs(g))<grad.tol)return(par)
    h <<- he(par)
    if(smartsearch){
      fnpar <- fn(par)
      p <- NULL
      f <- function(t,gradient=FALSE){
        ## Fast check: negative diagonal elements
        m <- min(diag(h))
        if(m<0){
          if(!(t>-m)){ ## h+t*I negative definite
            ustep <<- min(ustep,invphi(-m))
            return(NaN)
          }
        }
        ## Passed...
        ## Now do more expensive check...
        ##ok <- !is.character(try( .Call("destructive_CHM_update",L,h,as.double(t),PACKAGE="Matrix") , silent=silent))
        ok <- !is.character(try( updateCholesky(L,h,t) , silent=silent))
        if(!ok)return(NaN)
        dp <- as.vector(solve(L,g))
        p <<- par-dp
        ans <- fn(p)
        if(gradient)attr(ans,"gradient") <- sum(solve(L,dp)*gr(p))
        ans
      }
      
      ## Adaptive stepsize algorithm (smartsearch)
      phi <- function(u)1/u-1
      invphi <- function(x)1/(x+1)
      fu <- function(u){f(phi(u))}
      ## ========== Functions controling the algorithm
      ## Important requirements:
      ## 1. increase(u) and decrease(u) takes values in [0,1]
      ## 2. increase(u)>u and decrease(u)<u
      ## 3. increase(u)->1 when u->1
      ## 4. decrease(u)->0 when u->0
      ## Properties of algorithm:
      ## * ustep must converge towards 1 (because 1 <==> Positive definite hessian)

      ## power<1 - controls the boundary *repulsion*
      increase <- function(u)u0+(1-u0)*u^power
      ##decrease <- function(u)1-increase(1-u)
      ## Solve problem with accuracy when u apprach 0
      decrease <- function(u)ifelse(u>1e-10,1-increase(1-u),(1-u0)*power*u)
      ##plot(increase,0,1,ylim=c(0,1));plot(decrease,0,1,add=TRUE);abline(0,1)
      ustep <<- increase(ustep)
      repeat{
        fu.value <- fu(ustep)
        if(is.finite(fu.value)){
          eps <- sqrt(.Machine$double.eps)
          if(fu.value>fnpar+eps){
            if(ustep<=0)break; ## Avoid trap
            ustep <<- decrease(ustep)
          }
          else break;
        } else {
          if(ustep<=0)break; ## Avoid trap
          ustep <<- decrease(ustep)
        }
      }
      if(trace>=1)cat("value:", fu.value,"mgc:",max(abs(g)), "ustep:", ustep ,"\n")
      return(p)
    }
    dpar <- chol.solve(h,g) ## ordinary newton
    if(trace>=1)cat("mgc:",max(abs(g)) ,"\n")
    par - alpha * dpar
  }
  norm <- function(x){
    res <- sqrt(sum(x^2))
    res
  }
  fn.history <- numeric(maxit)
  fail <- 0
  for (i in seq(length=maxit)){
    parold <- par
    if(trace>=1)cat("iter:",i," ")
    par <- iterate(par)
    fn.history[i] <- fn(par)
    if(i>10){
      tail10 <- tail(fn.history[1:i],10)
      improve10 <- tail10[1] - tail10[length(tail10)]  
      if(improve10<tol10){
        if(trace>=1)cat("Not improving much - will try early exit...")
        pd <- iterate(par,pd.check=TRUE)
        if(trace>=1)cat("PD hess?:",pd,"\n")
        if(pd)break;
        fail <- fail+1
      }
    }
    if(norm(par-parold)<step.tol){
      break
    }
    if(fail>5){
      stop("Newton drop out: Too many failed attempts.")
    }
  }
  pd <- iterate(par,pd.check=TRUE)
  if(!pd)stop("Newton failed to find minimum.")
  names(par) <- nam
  value <- fn(par)
  g <- gr(par)
  if(trace>=1)cat("mgc:",max(abs(g)),"\n")
  list(par=par,value=value,gradient=g,hessian=h,iterations=i)
}


sparseHessianFun <- function(obj,skipFixedEffects=FALSE){
  r <- obj$env$random
  if(skipFixedEffects){
    ## Assuming that random effects comes first in parameter list, we can set
    ## skip <- as.integer(length(obj$env$par)-length(r)) ## ==number of fixed effects
    skip <- seq.int(length.out=length(obj$env$par))[-r]
  } else {
    ##skip <- as.integer(0)
    skip <- integer(0) ## <-- Empty integer vector
  }
  ## ptr.list
  ADHess <- .Call("MakeADHessObject2", obj$env$data, obj$env$parameters, 
                  obj$env$reportenv,
                  skip, ## <-- Skip this index vector of parameters
                  PACKAGE=obj$env$DLL
                  )
  ev <- function(par=obj$env$par).Call("EvalADFunObject", ADHess$ptr, par,
                   control = list(
                     order = as.integer(0),
                     hessiancols = integer(0),
                     hessianrows = integer(0),
                     sparsitypattern = as.integer(0),
                     rangecomponent = as.integer(1),
                     dumpstack=as.integer(0)),PACKAGE=obj$env$DLL)
  i=as.integer(attr(ADHess$ptr,"i"))
  j=as.integer(attr(ADHess$ptr,"j"))
  n <- length(obj$env$par)
  M <- new("dsTMatrix",i=i,j=j,x=ev(),Dim=as.integer(c(n,n)),uplo="L")
  Hfull <- as(M,"dsCMatrix")
  Hrandom <- Hfull[r,r,drop=FALSE]
  function(par=obj$env$par,random=FALSE){
    if(!random){
      Hfull@x[] <- ev(par)
      return(Hfull)
    } else {
      if(skipFixedEffects){
        return( .Call("setxslot",Hrandom,ev(par),PACKAGE="TMB") )
      }
      else {
        Hfull@x[] <- ev(par)
        return(Hfull[r,r])
      }
    }
  }
}

## Debugging utility: Check sparse hessian.
## By comparing with gradient differentiated in random direction.
checkSparseHessian <- function(obj,par=obj$env$last.par,
                               w = rnorm(length(par)), ## random direction
                               plot=TRUE,...){
  r <- obj$env$random
  w[-r] <- 0
  res1 <- obj$env$f(par, order = 1, type = "ADGrad", rangeweight = w)[r]
  res2 <- (obj$env$spHess(par)%*%w)[r]
  res <- list(x=res1,y=res2)
  if(plot){
    plot(res,...)
    abline(0,1,col="red")
  }
  invisible(res)
}

runSymbolicAnalysis <- function(obj){
  ok <- .Call("have_tmb_symbolic",PACKAGE="TMB")
  if(!ok){
    cat("note: tmb_symbolic not installed\n")
    return(NULL)
  }
  h <- obj$env$spHess(random=TRUE)
  h@x[] <- 0
  diag(h) <- 1
  L <- .Call("tmb_symbolic",h,PACKAGE="TMB")
  obj$env$L.created.by.newton <- L
  NULL
}
