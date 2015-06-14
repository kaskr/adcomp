/** \file
\ingroup Examples
\brief Shows multivariate distributions

*/

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  PARAMETER(dummy_par)         

  // Load namespace which contains the multivariate distributions
  using namespace density;
	
  Type res;                    // Dummy variable to which results will be asigned
  
  // Example: MVNORM_t  -----------------------------------------------------------

  matrix<Type> Sigma(3,3);
  Sigma.fill(0.1);                   // Fill the whole matrix
  Sigma.diagonal() *= 10.0;          // Multiply diagonal by 10 to positive definite Sigma
  vector<Type> x0(3);                // Point of evaluation
  x0.fill(0.0);                      // Initialize x0 to be zero
  MVNORM_t<Type> N_0_Sigma(Sigma);   // N_0_Sigma is now a Distribution  
  res = N_0_Sigma(x0);               // Evaluates (neg. log) density at x
  res = MVNORM(Sigma)(x0);           // Shorthand form of the above two lines
  N_0_Sigma.cov();                   // Returns covariance matrix (Sigma in this case)
  REPORT(N_0_Sigma.cov());           // Report back to R
  REPORT(MVNORM_t<Type>(Sigma).cov()); // Should return the same matrix as line above
  
  // Example: UNSTRUCTURED_CORR_t  -----------------------------------------------------
    
  
   // First: read docs for UNSTRUCTURED_CORR_t
  vector<Type> Lx(6);				  // Free parameters to be estimated 
  Lx.fill(2.5);
  UNSTRUCTURED_CORR_t<Type> nll(Lx);
  vector<Type> x1(4);
  x1.fill(2.0);
  res = nll(x1	);
  REPORT(nll.cov());

  
  // Example: AR1_t  -----------------------------------------------------
  
  // Univariate case
  int n = 10;                   // Number of time steps
  vector<Type> x2(n);           // Evaluation point, i.e. the time series       
  x2.fill(1.0);
  Type phi = 0.5;				// Autocorrelation
  N01<Type> nllN01;             // Distribution to be used at each time step
  AR1_t<N01<Type> > nll2(phi,nllN01);  // Create AR(1) process with N(0,1) distribution
  REPORT(nll2(x2));				// Evaluate neg. log density
  AR1(phi)(x2);                 // Equivalent to line above
  
  // Scale variances
  vector<Type> sds2(n);         // Vector of standard deviations (SD)
  sds2.fill(2.0);               // Set SD's  
  REPORT(VECSCALE(nll2,sds2)(x2));   // Evaluate neg. log density
  REPORT(VECSCALE(AR1(phi),sds2)(x2)); // Equivalent to line above   
  
  // Multivariate AR(1) process where the innovation vector is correlated
  int p=2;                      // dim(x)
  array<Type> x(p,n);           // Evaluation point       
  x.fill(1.0);
  REPORT(x);					// Reported back to R so that we can see layout of x

  vector<Type> unconstrained_params(p*(p-1)/2);
  unconstrained_params.fill(0.01);	// Low correlaiton, but this is not directly the correlation  
  REPORT(AR1(phi,UNSTRUCTURED_CORR(unconstrained_params))(x)); // nll 

  // Find nll for univariate time series. 
  // Do not add to nll for system of time series due to intra series correlation.
  REPORT(AR1(phi)(x.transpose().col(0)));	// First time series
  REPORT(AR1(phi)(x.transpose().col(1)));	// Second time series
    
  return Type(0.0);
}