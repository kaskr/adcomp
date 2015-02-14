// Leukemia example from Lindgren et al 2011, JRSS-B
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
  PARAMETER(log_alpha);  
  PARAMETER_VECTOR(x);  

  Type tau = exp(log_tau);
  Type kappa = exp(log_kappa);
  Type alpha = exp(log_alpha);

  Type nll = 0.0;

  SparseMatrix<Type> Q = Q_spde(spde,kappa);

  nll = GMRF(Q)(x);                              // Negative log likelihood

  vector<Type> Xbeta = X*beta;  
  for(int i=0; i<time.size(); i++){    
    Type eta = Xbeta(i) + x(meshidxloc(i))/tau;
    Type lambda = exp(eta);
    Type t_alpha = pow(time(i),alpha);
    Type S = exp(-lambda*t_alpha);               // Survival function
    Type f = lambda*alpha*t_alpha/time(i)*S;     // Densities

    // Likelihood contribution depends on truncation status
    if(notcens(i))
      nll -= log(f);
    else
      nll -= log(S); 
  }
    
  return nll;
}
