// Demonstrate user specified atomic functions (parallel version).
#include <TMB.hpp>

template<class Type>
vector<Type> dowork(vector<Type> x){
  int n=400;
  vector<Type> y(n);
  for(int i=0;i<n;i++)y[i]=dnorm(Type(i)/Type(n),x[0],x[1],0);
  Type s=0;
  for(int i=0;i<n;i++)s+=y[i];
  for(int i=0;i<n;i++)y[i]/=s;
  return y;
}
REGISTER_ATOMIC(dowork);


template<class Type>
Type objective_function<Type>::operator() ()
{
  #ifdef _OPENMP
  this->max_parallel_regions=omp_get_max_threads();
  #endif

  PARAMETER_ARRAY(x);
  int m=x.cols();
  Type res=0;
  int n=400;
  vector<Type> tmp(n);
  tmp.setZero();
  for(int i=0;i<m;i++){
    PARALLEL_REGION tmp+=dowork(vector<Type>(x.col(i)));
  }
  res=tmp.sum();
  return res;
}
