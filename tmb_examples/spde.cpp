// Illustration SPDE/INLA approach to spatial modelling via Matern correlation function
// Leukemia example from Lindgren et al 2011, JRSS-B
// http://www.r-inla.org/examples/case-studies/lindgren-rue-and-lindstrom-rss-paper-2011

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  DATA_VECTOR(time);
  DATA_IVECTOR(notcens);
  DATA_IVECTOR(meshidxloc);
  DATA_MATRIX(X);    
  
  DATA_STRUCT(spde,spde_t);
  
  PARAMETER_VECTOR(beta);      
  PARAMETER(log_tau);
  PARAMETER(log_kappa);
  PARAMETER(log_omega);  
  PARAMETER_VECTOR(x);  

  Type tau = exp(log_tau);
  Type kappa = exp(log_kappa);
  Type omega = exp(log_omega);  // Parameter of Weibull distribution

  Type nll = 0.0;

  SparseMatrix<Type> Q = Q_spde(spde,kappa);

  nll = GMRF(Q)(x);                              // Negative log likelihood

  vector<Type> Xbeta = X*beta;  
  for(int i=0; i<time.size(); i++){    
    Type eta = Xbeta(i) + x(meshidxloc(i))/tau;
    Type lambda = exp(eta);
    Type t_omega = pow(time(i),omega);
    Type S = exp(-lambda*t_omega);               // Survival function
    Type f = lambda*omega*t_omega/time(i)*S;     // Weibull density

    // Likelihood contribution depends on truncation status
    if(notcens(i))
      nll -= log(f);
    else
      nll -= log(S); 
  }
 
  double nu = 1.0;            // nu = alpha-d/2 = 2-1 by eqn (2) in Lindgren 
  Type rho = sqrt(8*nu)/kappa;  // Distance at which correlation has dropped to 0.1 (p.  4 in Lindgren)
  ADREPORT(rho);
    
  return nll;
}
