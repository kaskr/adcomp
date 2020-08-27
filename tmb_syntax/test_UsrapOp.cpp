#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(A);
  PARAMETER_VECTOR(u);
  DATA_VECTOR(y);
  PARAMETER(sd);
  Type f = 0;
  f -= dnorm(u, Type(0), Type(1), true).sum();
  matrix<Type> um = u.matrix();
  vector<Type> ypred = atomic::matmul( A, um );
  f -= dnorm(y, ypred, sd, true).sum();
  return f;
}
