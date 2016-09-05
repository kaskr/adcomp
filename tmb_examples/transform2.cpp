// Beta distributed random effects using copulas.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);

  PARAMETER(phi);
  PARAMETER(shape1);
  PARAMETER(shape2);
  PARAMETER(sd);
  PARAMETER_VECTOR(u);

  Type res = 0;
  res += density::AR1(phi)(u);
  vector<Type> unif = pnorm(u, Type(0), Type(1));
  vector<Type> x = qbeta(unif, shape1, shape2);
  res -= dnorm(y, x, sd, true).sum();
  return res;
}
