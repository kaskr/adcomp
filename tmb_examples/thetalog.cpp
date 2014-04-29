// Theta logistic population model from Pedersen et al 2012, Ecol. Modelling.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* Data section */
  DATA_VECTOR(Y);
  /* Parameter section */
  PARAMETER_VECTOR(X);
  PARAMETER(logr0);
  PARAMETER(logtheta);
  PARAMETER(logK);
  PARAMETER(logQ);
  PARAMETER(logR);
  /* Procedure section */
  Type r0=exp(logr0);
  Type theta=exp(logtheta);
  Type K=exp(logK);
  Type Q=exp(logQ);
  Type R=exp(logR);
  int timeSteps=Y.size();
  Type ans=0;
  for(int i=1;i<timeSteps;i++){
    Type m=X[i-1]+r0*(1.0-pow(exp(X[i-1])/K,theta));
    ans-=dnorm(X[i],m,sqrt(Q),true);
  }
  for(int i=0;i<timeSteps;i++){
    ans-=dnorm(Y[i],X[i],sqrt(R),true);
  }
  return ans;
}
