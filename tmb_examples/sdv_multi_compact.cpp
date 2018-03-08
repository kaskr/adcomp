// Compact version of sdv_multi
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

  using namespace density;

  Type g=0;
  vector<Type> sigma=exp(log_sigma);
  array<Type> ht=h.transpose(); // For row access
  vector<Type> sigma_init=sigma/sqrt(Type(1.0)-phi*phi); // Initial sd of AR1
  for(int j=0;j<p;j++)g += SCALE(AR1(phi(j)),sigma_init(j))(ht.col(j));

  // Likelihood contribution: observations
  UNSTRUCTURED_CORR_t<Type> neg_log_density(off_diag_x);
  for(int i=0;i<n;i++)
  {
    vector<Type> sigma_y = exp( Type(0.5) * (mu_x + h.col(i).vec() ));
    g += VECSCALE(neg_log_density,sigma_y)( y.col(i).vec() );
  }

  return g;
}
