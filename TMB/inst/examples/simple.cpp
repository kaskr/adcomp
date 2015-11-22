// Normal linear mixed model specified through sparse design matrices.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);         // Observations
  DATA_SPARSE_MATRIX(B);  // Random effect design matrix
  DATA_SPARSE_MATRIX(A);  // Fixed effect design matrix
  PARAMETER_VECTOR(u);    // Random effects vector
  PARAMETER_VECTOR(beta); // Fixed effects vector
  PARAMETER(logsdu);      // Random effect standard deviations
  PARAMETER(logsd0);      // Measurement standard deviation

  // Distribution of random effect (u):
  Type ans = 0;
  ans -= dnorm(u, Type(0), exp(logsdu), true).sum();

  // Distribution of obs given random effects (x|u):
  vector<Type> y = A * beta + B * u;
  ans -= dnorm(x, y, exp(logsd0), true).sum();

  // Apply delta method on sd0:
  ADREPORT( exp(logsd0) );

  // Report posterior mode and mean of sum(exp(u))
  ADREPORT( sum(exp(u)) );

  return ans;
}

