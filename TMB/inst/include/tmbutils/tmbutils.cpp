/** \file 
   \brief Namespace of utility functions for TMB
*/
namespace tmbutils{
// Utilities used by the core
using namespace Eigen;
#include "vector.cpp"
#include "array.cpp"

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

