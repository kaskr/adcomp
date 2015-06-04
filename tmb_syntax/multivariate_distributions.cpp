// Shows use of multivariate distributions
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
  int p=2;                      // dim(x)
  array<Type> x(p,n);           // Evaluation point       
  x.fill(1.0);
  REPORT(x);					// Reported back to R so that we can see layout of x

  vector<Type> unconstrained_params(p*(p-1)/2);
  unconstrained_params.fill(0.01);	// Low correlaiton, but this is not directly the correlation  
  Type phi = 0.5;				// Autocorrelation
  REPORT(AR1(phi,UNSTRUCTURED_CORR(unconstrained_params))(x)); // nll 

  // Find nll for univariate time series. 
  // Do not add to nll for system of time series due to intra series correlation.
  REPORT(AR1(phi)(x.transpose().col(0)));	// First time series
  REPORT(AR1(phi)(x.transpose().col(1)));	// Second time series
    
  return Type(0.0);
}

/** \file
\ingroup Examples
\brief Shows multivariate distributions

*/
