/* Templates to get convenient R-like syntax. */


/* General utility similar to R's split function. */
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

/* Sum of vector, matrix or array */
template<template<class> class Vector, class Type>
Type sum(Vector<Type> x){return x.sum();}
