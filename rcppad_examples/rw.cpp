// Random walk with multivariate correlated increments and measurement noise.
#include <RcppAD.hpp>

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
  int timeSteps=obs.dim[0];
  int stateDim=obs.dim[1];
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
  ans-=dnorm(vector<Type>(u(0)),Type(0),Type(1),1).sum();
  for(int i=1;i<timeSteps;i++)    
    ans+=neg_log_density(u(i)-u(i-1)); // Process likelihood 
  for(int i=1;i<timeSteps;i++)
    ans-=dnorm(vector<Type>(obs(i)),vector<Type>(u(i)),sdObs,1).sum(); // Data likelihood
  return ans;
}
