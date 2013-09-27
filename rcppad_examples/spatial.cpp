// Spatial poisson GLMM on a grid, with exponentially decaying correlation function
#include <RcppAD.hpp>

/* Parameter transform */
template <class Type>
Type Exp(Type x){return exp( x);}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n);
  DATA_VECTOR(y);
  DATA_MATRIX(X)
  DATA_MATRIX(dd)
  PARAMETER_VECTOR(b);
  PARAMETER(a);
  PARAMETER(log_sigma);
  PARAMETER_VECTOR(u);
  Type sigma2=Exp(2.0*log_sigma);

  using namespace density;
  int i,j;
  Type res=0;

  vector<Type> eta(n); 
  eta  =   X*b + u;

  // 
  matrix<Type> cov(n,n); 
  for (i=0;i<n;i++)
  {
    cov(i,i)=sigma2;
    for ( j=0;j<i;j++)
    {
      cov(i,j)=sigma2*exp(-a*dd(i,j));			// Exponentially decaying correlation
      cov(j,i)=cov(i,j);
    }
  }

  MVNORM_t<Type> neg_log_density(cov);
  res+=neg_log_density(u);
  

  // logdpois = N log lam - lam
  for(i=0;i<n;i++) res -= y[i]*eta[i]-exp(eta[i]);

  return res;

}
