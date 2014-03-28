// socatt from ADMB example collection.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  DATA_FACTOR(y); //categorical response vector
  DATA_INTEGER(S); //number of response categories
  DATA_MATRIX(X); // Fixed effects design matrix
  DATA_FACTOR(group);
  PARAMETER_VECTOR(b); // Fixed effects
  PARAMETER(logsigma);
  PARAMETER_VECTOR(tmpk); // kappa ( category thresholds)
  PARAMETER_VECTOR(u);    // Random effects
  Type sigma = exp(logsigma);
  vector<Type> alpha = tmpk;
  for(int s=1;s<tmpk.size();s++)
    alpha(s) = alpha(s-1) + exp(tmpk(s));
  Type ans=0;
  ans -= sum(dnorm(u,Type(0),Type(1),true));
  vector<Type> eta = X*b;
  for(int i=0; i<y.size(); i++){
    eta(i) += sigma*u(group(i));
    Type P;
    if(y(i)==(S-1)) P = 1.0; else P = Type(1)/(Type(1)+exp(-(alpha(y(i))-eta(i))));
    if(y(i)>0) P -= Type(1)/(Type(1)+exp(-(alpha(y(i)-1)-eta(i))));
    ans -= log(1.e-20+P);
  }
  
  return ans;
}
