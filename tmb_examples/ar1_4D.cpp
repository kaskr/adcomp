// Separable covariance on 4D lattice with AR1 structure in each direction.
#include <TMB.hpp>

/* Parameter transform */
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(N)
  PARAMETER_ARRAY(eta);
  PARAMETER(transf_phi); /* fastest running dim */
  Type phi=f(transf_phi);
  ADREPORT(phi);

  using namespace density;
  Type res=0;

  res+=AR1(phi,AR1(phi,AR1(phi,AR1(phi))))(eta);

  // logdpois = N log lam - lam
  for(int i=0;i<N.size();i++)res-=N[i]*eta[i]-exp(eta[i]);

  SIMULATE {
    AR1(phi,AR1(phi,AR1(phi,AR1(phi)))).simulate(eta);
    vector<Type> lam = exp(eta);
    N = rpois(lam);
    REPORT(eta);
    REPORT(N);
  }

  return res;

}
