// Conway-Maxwell-Poisson distribution
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);
  PARAMETER(logmu); Type mu = exp(logmu);
  PARAMETER(lognu); Type nu = exp(lognu);
  DATA_STRING(parameterization);

  Type f = 0;

  if (parameterization == "mode") {
    Type lambda = exp(nu * log(mu));
    f -= sum(dcompois(x, lambda, nu, true));
  }
  else
  if (parameterization == "mean") {
    f -= sum(dcompois2(x, mu, nu, true));
  }
  else
    error("Unknown parameterization");

  ADREPORT(mu);
  ADREPORT(nu);

  return f;
}
