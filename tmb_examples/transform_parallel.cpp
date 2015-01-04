// Parallel version of transform
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);

  PARAMETER(phi);
  PARAMETER(shape);
  PARAMETER(scale);
  PARAMETER(sd);
  PARAMETER_VECTOR(u);

  parallel_accumulator<Type> res(this);
  res += density::AR1(phi)(u);
  vector<Type> unif = pnorm(u,Type(0),Type(1));
  vector<Type> x = qgamma(unif,shape,scale);
  for(int i=0;i<x.size();i++)res -= dnorm(y[i],x[i],sd,true);
  return res;
}

