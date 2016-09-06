// Gaussian process with Matern covariance.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);
  DATA_MATRIX(D);
  PARAMETER(phi);
  PARAMETER(kappa);

  matrix<Type> C(D);
  for(int i=0; i<C.rows(); i++)
    for(int j=0; j<C.cols(); j++)
      C(i,j) = matern(D(i,j), phi, kappa);

  Type nll = density::MVNORM_t<Type>(C)(x);

  return nll;
}
