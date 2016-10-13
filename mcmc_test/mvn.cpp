#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(d); // dimension of matrix
  DATA_MATRIX(covar);
  PARAMETER_VECTOR(X); // parameters
  using namespace density;
  MVNORM_t<Type> neg_log_density(covar);
  Type nll= neg_log_density(X);
  REPORT(nll);
  return(nll);
}
