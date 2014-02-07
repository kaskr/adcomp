// Random walk with multivariate correlated increments and measurement noise.
#include <TMB.hpp>
 
/* Parameter transform */
template <class Type>
Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y); /* timeSteps x stateDim */
  PARAMETER_VECTOR(X); /* State */
  PARAMETER(logr0);
  PARAMETER(logtheta);
  PARAMETER(logK);
  PARAMETER(logQ)
  PARAMETER(logR)
  
  int timeSteps=Y.size();
  using namespace density;
  Type ans=0;
  for(int i=1;i<timeSteps;i++){
    Type var=exp(logQ);
    Type m=X[i-1] + exp(logr0) * (1.0 - pow(exp(X[i-1])/exp(logK),exp(logtheta))); 
    ans+=0.5*(log(2.0*M_PI*var)+square(X[i]-m)/var);
  }
  for(int i=0;i<timeSteps;i++){
    Type var=exp(logR);
    ans+=0.5*(log(2.0*M_PI*var)+square(X[i]-Y[i])/var);
  }
  return ans;
}
