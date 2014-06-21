// Collection of .Call()'s in "TMB.R". The purpose of this file is to show which C++ functions are called from R.
// Each function call is proceeded by the corresponding line from TMB.R
// Warning: this file may be outdated. 
#include <TMB.hpp>

//  .Call("destructive_CHM_update",L,H,as.double(t),PACKAGE="Matrix")
destructive_CHM_update();

//    parNameOrder <- .Call("getParameterOrder",data,parameters,new.env(),PACKAGE=DLL)
getParameterOrder();

//      ADFun <<- .Call("MakeADFunObject",data,parameters,reportenv,
MakeADFunObject();

//      Fun <<- .Call("MakeDoubleFunObject",data,parameters,reportenv,PACKAGE=DLL)
MakeDoubleFunObject();

//      .Call("EvalDoubleFunObject",Fun$ptr,unlist(parameters),control=list(order=as.integer(0)),PACKAGE=DLL)
EvalDoubleFunObject();

//      ADGrad <<- .Call("MakeADGradObject",data,parameters,reportenv,PACKAGE=DLL)
MakeADGradObject();

//      res <- .Call("EvalADFunObject",ADFun$ptr,theta,
//        .Call("EvalADFunObject",e$ADHess$ptr,theta,
//  ev <- function(par=obj$env$par).Call("EvalADFunObject", ADHess$ptr, par,
EvalADFunObject();

//      res <- .Call("EvalDoubleFunObject",Fun$ptr,theta,
EvalDoubleFunObject();

//      solveSubset <- function(L).Call("tmb_invQ",L,PACKAGE="TMB")
tmb_invQ();

//      solveSubset2 <- function(L).Call("tmb_invQ_tril_halfdiag",L,PACKAGE="TMB")
tmb_invQ_tril_halfdiag();

//        m <- .Call("match_pattern",A,B,PACKAGE="TMB") ## Same length as A@x with pointers to B@x
match_pattern();

//      ##.Call("destructive_CHM_update",L,hessian,as.double(0),PACKAGE="Matrix")
destructive_CHM_update();

//  .Call("omp_num_threads",n,PACKAGE="TMB")
omp_num_threads();

//    unlist(.Call("InfoADFunObject",get(name,env),PACKAGE=obj$env$DLL))
InfoADFunObject();

//    unlist(.Call("optimizeADFunObject",get(name,env)$ptr,PACKAGE=obj$env$DLL))
//    .Call("optimizeADFunObject",ADHess$ptr,PACKAGE=obj$env$DLL)
optimizeADFunObject();

//  ADHess <- .Call("MakeADHessObject2", obj$env$data, obj$env$parameters, 
MakeADHessObject2();

//        return( .Call("setxslot",Hrandom,ev(par),PACKAGE="TMB") )
setxslot();

//  ok <- .Call("have_tmb_symbolic",PACKAGE="TMB")
have_tmb_symbolic();

//  L <- .Call("tmb_symbolic",h,PACKAGE="TMB")
tmb_symbolic();
