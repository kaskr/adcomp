// Shows use of multivariate distributions
// Objects of "double" type (instead of Type) so that they can be printed
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Input from R (not used)
  PARAMETER(dummy_par)

  // Load namespace which contains the multivariate distributions
  using namespace density;
  
  // Multivariate AR(1) process where the innovation vector is correlated
  int n = 10;                   // Number of time steps
  int p=3;                      // dim(x)
  array<Type> x(n,p);           // Evaluation point       
  x.fill(1.0);
  REPORT(x);
  vector<double> unconstrained_params(p*(p-1)/2);
  unconstrained_params.fill(0.01);	// This is not the correlation itself!!
  double phi = 0.5;				// Autocorrelation
  //AR1(phi,UNSTRUCTURED_CORR(unconstrained_params))(x);
  

  Type ans = Type(0.0);
  return ans;

}

/** \file
\ingroup Examples
\brief Shows multivariate distributions

*/
