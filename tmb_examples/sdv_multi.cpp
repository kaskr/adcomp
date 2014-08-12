// Multivatiate SV model from Skaug and Yu 2013, Comp. Stat & data Analysis (to appear)
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n);
  DATA_INTEGER(p);
  DATA_ARRAY(y);
  PARAMETER_VECTOR(phi);
  PARAMETER_VECTOR(log_sigma);
  PARAMETER_VECTOR(mu_x);
  PARAMETER_VECTOR(off_diag_x);
  PARAMETER_ARRAY(h);

  int i,j;
  Type g=0;

  vector<Type> sigma=exp(log_sigma);

  // Likelihood contribution: stationary distribution for initial state
  vector<Type> tmp(p);
  for(j=0;j<p;j++)
    tmp(j) = sigma(j)/sqrt(Type(1.0)-phi(j)*phi(j)); 
  g -= sum(dnorm(vector<Type>(h.col(0)),Type(0),tmp,1));

  // Likelihood contribution: State transitions
  for(i=1;i<n;i++)
  {
    vector<Type> tmp2(p);
    for(j=0;j<p;j++)
      tmp2(j) = phi(j)*h(j,i-1);
    g -= sum(dnorm(vector<Type>(h.col(i)),tmp2,sigma,1));
  }

  // Cholesky factor of Sigma
  matrix<Type> L(p,p);
  L.setIdentity();

  int k=0;
  for(i=1;i<p;i++)
  {
    Type Norm2=L(i,i);
    for(j=0;j<=i-1;j++)
    {
      L(i,j) = off_diag_x(k++);
      Norm2 += L(i,j)*L(i,j); 
    }
    for(j=0;j<=i;j++)
      L(i,j) /= sqrt(Norm2);
  }

  matrix<Type> Sigma = L * L.transpose();

  using namespace density;

  // Likelihood contribution: observations
  for(i=0;i<n;i++)
  {
    vector<Type> sigma_y = exp(Type(0.5)*(mu_x + vector<Type>(h.col(i))));

    // Scale up correlation matrix
    matrix<Type> Sigma_y(p,p);
    for(int i2=0;i2<p;i2++)
      for(j=0;j<p;j++)
        Sigma_y(i2,j) = Sigma(i2,j)*sigma_y(i2)*sigma_y(j);

    g += MVNORM(Sigma_y)(vector<Type>(y.col(i)));
  }

  return g;
}
