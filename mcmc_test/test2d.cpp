// Test distributional properties of 2D MCMC samples
#include <TMB.hpp>

/* 
   Specify a model through factors of joint 2D CDF:
   P(X <= x) and 
   P(Y <= y | X = x)

   This allow us to

   1. Use autodiff to get the joint density, which is passed to MCMC
   sampler.

   2. Easily transform MCMC samples to IID uniform variates.
*/
struct JointCDF {
  int example;
  template <class T>
  vector<T> operator()(vector<T> xy){
    T x = xy(0), y = xy(1);
    vector<T> ans(2);
    switch(example) {
    case 0 :
      ans <<
	pnorm( x + sin( 0.8 * x) )  ,
	pnorm( y - sin( 3.0 * x) )   ;
      break;
    default:
      error("Not implemented");
    }
    return ans;
  }
};

template<class Type>
Type objective_function<Type>::operator() ()
{
  PARAMETER(x);
  PARAMETER(y);
  DATA_INTEGER(example);
  vector<Type> xy(2);
  vector<Type> uv(2);
  xy << x, y;
  JointCDF F = { example };
  matrix<Type> J = autodiff::jacobian(F, xy);
  uv = F(xy);
  Type f = -log(J(0,0)) - log(J(1,1));
  REPORT(uv);
  return f;
}
