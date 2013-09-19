//  --------------------------------------------------------------------------
// Copyright (c) 2008,2009,2010,2011,2012, Anders Nielsen <an@aqua.dtu.dk> 
// and Casper Berg <cbe@aqua.dtu.dk>. All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//   * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//   * Neither the name of the assessment tool SAM nor the
//     names of its contributors may be used to endorse or promote products
//     derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
// ARE DISCLAIMED. IN NO EVENT SHALL ANDERS NIELSEN OR CASPER BERG BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY 
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH 
// DAMAGE.
//  --------------------------------------------------------------------------

#ifndef __nLogNormal_h__
#define __nLogNormal_h__

#include <fvar.hpp>
#include <math.h>

_CONST double log2pi = log(2.0*M_PI);

template<class MatrixType>
MatrixType bksolve(_CONST MatrixType& L, _CONST MatrixType& b)
{
  RETURN_ARRAYS_INCREMENT();
  int i, k, r, R, c, C;
  r=b.rowmin();
  R=b.rowmax();
  c=b.colmin();
  C=b.colmax();
  MatrixType sumVec(1,1,c,C);
  MatrixType x(r,R,c,C);
  for(i=r; i<=R; ++i){
    sumVec(1)=b(i);
    for(k=i-1; k>=r; k--){
      sumVec(1)-=L(i,k)*x(k);
    }
    x(i)=sumVec(1)/L(i,i);
  }
  RETURN_ARRAYS_DECREMENT();
  return x;
}

template <class VectorType, class MatrixType>
  VectorType nLogNormal(_CONST VectorType& x, _CONST MatrixType& mu, _CONST MatrixType& S)
{
  RETURN_ARRAYS_INCREMENT();
  int r, R, c, C, N;
  r=mu.rowmin();
  R=mu.rowmax();
  c=mu.colmin();
  C=mu.colmax();
  N=R-r+1;
  VectorType ret(c,C);
  MatrixType diff(r,R,c,C);
  for(int i=c; i<=C; ++i){
    for(int j=r; j<=R; ++j){
      diff(j,i)=x(j)-mu(j,i);
    }
  }
  VectorType logDet(1,1);
  logDet(1)=0.0;
  MatrixType chol=choleski_decomp(S);
  for(int i=r; i<=R; ++i){logDet(1)+=log(chol(i,i));}
  logDet(1)*=2.0;
  MatrixType tmp=bksolve(chol,diff);
  ret=0.5*(log2pi*N+logDet(1)+colsum(square(tmp)));
  RETURN_ARRAYS_DECREMENT();
  return ret;
}

template <class VectorType, class MatrixType>
  VectorType nLogNormal(_CONST MatrixType& x, _CONST VectorType& mu, _CONST MatrixType& S)
{
  RETURN_ARRAYS_INCREMENT();
  VectorType ret=nLogNormal(mu, x, S);
  RETURN_ARRAYS_DECREMENT();
  return ret;
}

df1b2variable nLogNormal(_CONST df1b2vector& x, _CONST df1b2vector& mu, _CONST df1b2matrix& S)
{
  RETURN_ARRAYS_INCREMENT();
  int r=mu.indexmin(), R=mu.indexmax();
  df1b2matrix MU(r,R,1,1);
  for(int i=r; i<=R; ++i){MU(i,1)=mu(i);}
  df1b2vector tmp=nLogNormal(x, MU, S);
  df1b2variable ret=tmp(tmp.indexmin());
  RETURN_ARRAYS_DECREMENT();
  return ret;
}

dvariable nLogNormal(_CONST dvar_vector& x, _CONST dvar_vector& mu, _CONST dvar_matrix& S)
{
  RETURN_ARRAYS_INCREMENT();
  dvar_matrix MU(mu.indexmin(),mu.indexmax(),1,1);
  MU.colfill(1,mu);
  dvar_vector tmp=nLogNormal(x, MU, S);
  dvariable ret=tmp(tmp.indexmin());
  RETURN_ARRAYS_DECREMENT();
  return ret;
}

#endif



 
