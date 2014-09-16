/** \file 
   \brief Namespace of utility functions for TMB
*/
namespace tmbutils{
//#include <Eigen/Dense>
using namespace Eigen;
#include "vector.cpp"
#include "array.cpp"
#include "spmat.cpp"

#include "kronecker.cpp"
#include "matexp.cpp"
#include "splines.cpp"
#include "order.cpp"

template<class Type, class T1, class T2>
vector<Type> dnorm(vector<Type> x, T1 mean, T2 sd, int give_log=0)
{
  vector<Type> logres;
  x=(x-mean)/sd;
  logres=-log(Type(sqrt(2*M_PI))*sd)-Type(.5)*x*x;
  if(give_log)return logres; else return exp(logres);
}
  
template <class Type, class From>
vector<Type> asVector(From *px, int n){
  vector<Type> x(n);
  for(int i=0;i<n;i++)x[i]=Type(px[i]);
  return x;
}

#if defined(R_R_H) 
template <class Type>
array<Type> asArray(SEXP x)
{
  if(!isArray(x))error("NOT AN ARRAY!");
  SEXP dim=getAttrib(x,R_DimSymbol);
  vector<int> d=asVector<int,int>(INTEGER(dim), LENGTH(dim));
  vector<Type> y=asVector<Type,double>(REAL(x), LENGTH(x));
  return array<Type>(y,d);
}


#endif


} // End namespace

/** 
   \brief Collection of multivariate Gaussian distributions (members listed in \ref density.cpp)

   \ingroup Densities
   
   For use of the namespace see \ref Densities
 
   
*/
namespace density{
  using namespace tmbutils;
#include "density.cpp"
} // End namespace
