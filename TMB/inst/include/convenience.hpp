/** \file
  \brief Templates to get convenient R-like syntax. 
*/

/** \brief  Similar to R's split function:  split(x,fac) devides  x  into groups defined by  fac . 
* \details
Returns a "vector of vectors". 
*/
template<class Type>
vector<vector<Type> > split(vector<Type> x,vector<int> fac){
  if(x.size()!=fac.size())error("x and fac must have equal length.");
  int nlevels=0;
  for(int i=0;i<fac.size();i++)if(fac[i]>=nlevels)nlevels++;
  vector<vector<Type> > ans(nlevels);
  vector<int> lngt(nlevels);
  lngt.setZero();
  for(int i=0;i<fac.size();i++)lngt[fac[i]]++;
  for(int i=0;i<nlevels;i++)ans[i].resize(lngt[i]);
  lngt.setZero();
  for(int i=0;i<fac.size();i++){
    ans[fac[i]][lngt[fac[i]]]=x[i];
    lngt[fac[i]]++;
  }
  return ans;
}

/** Sum of vector, matrix or array */
template<template<class> class Vector, class Type>
Type sum(Vector<Type> x){return x.sum();}

/** Matrix * vector 

  Simplifies syntax in that .matrix() can be avoided. Recall: TMB type vector is of Eigen type Array.
*/
template<class Type>
vector<Type> operator*(matrix<Type> A, vector<Type> x){return A*x.matrix();}

/** SparseMatrix * vector */
template<class Type>
vector<Type> operator*(Eigen::SparseMatrix<Type> A, vector<Type> x){return (A*x.matrix()).array();}
