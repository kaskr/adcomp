#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(a);
  PARAMETER_VECTOR(x);

  using namespace atomic;

  Type res=0;
  if(a==0){
    for(int i=0;i<x.size();i++)res+=pnorm(x[i]);
  } else if(a==1){
    for(int i=0;i<x.size();i++)res+=qnorm(x[i]);
  } else if(a==2){
    CppAD::vector<Type> arg(3);
    arg[0] = x[0];
    arg[1] = x[1];
    arg[2] = Type(0);
    res += D_incpl_gamma_shape(arg)[0];
  }
  else if(a==3){
    CppAD::vector<Type> xx(x.size());
    for(int i=0;i<x.size();i++)xx[i]=x[i];
    vector<Type> y=matmul(xx);
    res = y.sum();
  }
  else if(a==4){
    CppAD::vector<Type> xx(x.size());
    for(int i=0;i<x.size();i++)xx[i]=x[i];
    vector<Type> y=matinv(xx);
    res = y.sum();
  }
  else if(a==5){
    CppAD::vector<Type> xx(x.size());
    for(int i=0;i<x.size();i++)xx[i]=x[i];
    res = logdet(xx)[0];
  }
  else if(a==6){
    CppAD::vector<Type> arg(2);
    arg[0] = x[0];
    arg[1] = x[1];
    res += inv_incpl_gamma(arg)[0];
  }
  else if(a==7){
    CppAD::vector<Type> arg(2);
    arg[0] = x[0];
    arg[1] = Type(0);
    res += D_lgamma(arg)[0];
  }
  else if(a==8){
    res += pgamma(x[0],x[1],x[2]);
  }
  else if(a==9){
    res += qgamma(x[0],x[1],x[2]);
  }
  else {
    error("Invalid a");
  }
  return res;
}

