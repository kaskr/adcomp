#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(obs);
  DATA_FACTOR(group);
  PARAMETER_VECTOR(mu);
  PARAMETER_VECTOR(sd);
  Type res=0;
  for(int i=0;i<obs.size();i++){
    res -= dnorm(obs[i],mu[group[i]],sd[group[i]],true);
  }
  return res;
}

