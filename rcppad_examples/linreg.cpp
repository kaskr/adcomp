#include <RcppAD.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  DATA_VECTOR(x);
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(logSigma);
  Type nll=-dnorm(Y,a+b*x,exp(logSigma),true).sum();
  return nll;
}

