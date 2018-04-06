// Identical with random walk example. Utilizing sparse block structure so efficient when the number of states is high.
#include <TMB.hpp>

/* Parameter transform */
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_ARRAY(obs); /* timeSteps x stateDim */
  PARAMETER_ARRAY(u); /* State */
  PARAMETER(transf_rho);
  PARAMETER_VECTOR(logsds);
  PARAMETER_VECTOR(logsdObs);
  int timeSteps=obs.dim[1];
  Type rho=f(transf_rho);
  vector<Type> sds=exp(logsds);
  vector<Type> sdObs=exp(logsdObs);
  // Setup object for evaluating multivariate normal likelihood
  using namespace density;
  VECSCALE_t<AR1_t<N01<Type> > > neg_log_density=VECSCALE(AR1(rho),sds);
  /* Define likelihood */
  Type ans=0;
  ans -= dnorm(u.col(0).vec(), Type(0), Type(1), true).sum();
  for(int i=1; i<timeSteps; i++)
    ans += neg_log_density(u.col(i)-u.col(i-1)); // Process likelihood
  for(int i=1; i<timeSteps; i++)
    ans -= dnorm(obs.col(i).vec(), u.col(i).vec(), sdObs, true).sum(); // Data likelihood
  return ans;
}
