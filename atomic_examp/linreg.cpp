#include <TMB.hpp>

template<class Type>
vector<Type> hej(vector<Type> x){
  vector<Type> y(1);
  y(0) = dnorm(x, Type(0), Type(1), true).sum();
  return y;
}
REGISTER_ATOMIC(hej)

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  DATA_VECTOR(x);
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(logSigma);
  ADREPORT(exp(2*logSigma));
  //Type nll=-sum(dnorm(Y,a+b*x,exp(logSigma),true));
  vector<Type> tmp = (Y - (a+b*x));
  Type nll = hej(tmp)[0] +  logSigma*logSigma;
  return nll;
}

