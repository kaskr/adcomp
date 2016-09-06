// Scaled up version of the Orange Tree example (5000 latent random variables)
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_VECTOR(t);
  DATA_INTEGER(M);
  DATA_FACTOR(ngroup);
  DATA_INTEGER(multiply);
  PARAMETER_VECTOR(beta);
  PARAMETER(log_sigma);
  PARAMETER(log_sigma_u);
  PARAMETER_VECTOR(u);

  Type sigma=exp(log_sigma);
  Type sigma_u=exp(log_sigma_u);

  ADREPORT(sigma);
  ADREPORT(sigma_u);

  using namespace density;
  int i,j,k,ii;

  Type g=0;

  for(k=0;k< multiply;k++)
  {
    ii = 0;
    for(i=0;i< M;i++)
    {
      // Random effects contribution
      Type u1 = u[i+k*M];
      g -= -(log_sigma_u);
      g -= -.5*pow(u1/sigma_u,2);

      vector<Type> a(3); 
      a[0] = 192.0 + beta[0] + u1;
      a[1] = 726.0 + beta[1];
      a[2] = 356.0 + beta[2];

      Type tmp;
      Type f;

      for(j=0;j<ngroup(i);j++)
      {
        f = a[0]/(1+exp(-(t[ii]-a[1])/a[2]));
        tmp = (y[ii] - f)/sigma;
        g -= -log_sigma - 0.5*tmp*tmp;
        ii++;
      }
    }
  }

  return g;
}
