#include "TMB.h"
#include "distrib.h"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);
  DATA_VECTOR(n);
  DATA_MATRIX(A);
  PARAMETER_VECTOR(b);
  vector<Type> mu = A * b;
  PARAMETER(logsd);
  Type sd = exp(logsd);
  Type ans = 0;
  Type tiny = 0.0; // Set to 1e-12 for robustness
  for(int i=0; i < x.size(); i++) {
    ans -= log( my_atomic::GaussBinomial(x(i), n(i), mu(i), sd) + tiny );
  }
  return ans;
}
