\name{sdreport}
\alias{sdreport}
\title{General sdreport function.}
\usage{
  sdreport(obj, par.fixed = NULL, hessian.fixed = NULL)
}
\arguments{
  \item{obj}{Object returned by \code{MakeADFun}}

  \item{par.fixed}{Optional. Fixed effect parameter
  estimate (will be known to \code{obj} when an
  optimization has been carried out).}

  \item{hessian.fixed}{Optional. Hessian wrt. fixed effects
  (will be calculated from \code{obj} if missing).}
}
\value{
  Object of class \code{sdreport}
}
\description{
  After optimization of an AD model, \code{sdreport} is
  used to calculate standard deviations of all model
  parameters, including non linear functions of random and
  fixed effects parameters specified through the ADREPORT()
  macro from the user template.
}
\details{
  First, the Hessian wrt. the fixed effect parameter vector
  (\eqn{\theta}) is calculated. The fixed effects
  covariance matrix is approximated by
  \deqn{V(\hat\theta)=\nabla^2 l(\hat\theta)^{-1}} where
  \eqn{l} denotes the likelihood function (i.e.
  \code{obj$fn}).

  For non-random effect models the standard delta-method is
  used to calculate the covariance matrix. Let
  \eqn{\phi(\theta)} denote some non-linear function of
  \eqn{\theta}. Then \deqn{V(\phi(\hat\theta))\approx
  \nabla\phi V(\hat\theta) \nabla\phi'}

  For random effect models a generalized delta-method is
  used. First the joint covariance of random and fixed
  effects is estimated by \deqn{V \pmatrix{ \hat u \cr
  \hat\theta } \approx \pmatrix{ H_{uu}^{-1} & 0 \cr 0 & 0
  } + J V(\hat\theta) J' } where \eqn{H_{uu}} denotes
  random effect block of the full joint Hessian and \eqn{J}
  denotes the Jacobian of \eqn{\pmatrix{\hat u(\theta) \cr
  \theta}} wrt. \eqn{\theta}. Here, the first term
  represents the expected conditional variance given the
  fixed effects and the second term represents the variance
  of the conditional mean wrt. the fixed effects. Now the
  delta method can be applied on a general non-linear
  function \eqn{\phi(u,\theta)} of random effects \eqn{u}
  and fixed effects \eqn{\theta}: \deqn{V(\phi(\hat
  u,\hat\theta))\approx \nabla\phi V \pmatrix{ \hat u \cr
  \hat\theta }\nabla\phi'}
}

