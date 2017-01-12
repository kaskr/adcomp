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
    f -= sum(dcompois(x, mu /* mode */, nu, true));
  }
  else
  if (parameterization == "mean") {
    f -= sum(dcompois2(x, mu /* mean */, nu, true));
  }
  else
    error("Unknown parameterization");

  ADREPORT(mu);
  ADREPORT(nu);

  return f;
}
