// Random walk with multivariate correlated increments and measurement noise.
#include <TMB.hpp>

/* Parameter transform */
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_ARRAY(obs); /* timeSteps x stateDim */
  PARAMETER(transf_rho);
  PARAMETER_VECTOR(logsds);
  PARAMETER_VECTOR(logsdObs);
  PARAMETER_ARRAY(u); /* State */
  int timeSteps=obs.dim[1];
  int stateDim=obs.dim[0];
  Type rho=f(transf_rho);
  vector<Type> sds=exp(logsds);
  vector<Type> sdObs=exp(logsdObs);
  // Setup object for evaluating multivariate normal likelihood
  matrix<Type> cov(stateDim,stateDim);
  for(int i=0;i<stateDim;i++)
    for(int j=0;j<stateDim;j++)
      cov(i,j)=pow(rho,Type(abs(i-j)))*sds[i]*sds[j];
  using namespace density;
  MVNORM_t<Type> neg_log_density(cov);
  /* Define likelihood */
  Type ans=0;
  for(int i=1;i<timeSteps;i++)    
    ans += neg_log_density(u.col(i)-u.col(i-1)); // Process likelihood
  for(int i=0; i<timeSteps; i++)
    ans -= dnorm(obs.col(i).vec(), u.col(i).vec(), sdObs, true).sum(); // Data likelihood
  return ans;
}
