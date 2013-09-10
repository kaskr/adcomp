#include <RcppAD.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  PARAMETER_VECTOR(x);
  Type res=0;
  for(int i=0;i<x.size();i++)res+=x[i];
  res=res*res;
  return res;
}
