// Separable covariance on lattice with AR1 structure in each direction.
#include <RcppAD.hpp>

/* Parameter transform */
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(N)
  PARAMETER_ARRAY(eta);
  PARAMETER(transf_phi1); /* fastest running dim */
  PARAMETER(transf_phi2); /* slowest running dim */
  Type phi1=f(transf_phi1);
  Type phi2=f(transf_phi2);

  using namespace density;
  Type res=0;
  // phi1 fastest running
  // res+=AR1(phi2,AR1(phi1))(eta);
  // Equivalent:
  res+=SEPARABLE(AR1(phi2),AR1(phi1))(eta);

  // logdpois = N log lam - lam
  for(int i=0;i<N.size();i++)res-=N[i]*eta[i]-exp(eta[i]);

  return res;

}
