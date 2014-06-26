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


/** \brief  Approximate normal cumulative distribution function, similar to R's pnorm (one-argument case only).
* \details
To be replaced by more accurate version based on Rmath library.
*/
template<class Type>
Type pnorm_approx(Type x){
  Type a=993./880.;
  Type b=89./880.;
  x = x/sqrt(Type(2));
  return Type(.5) * tanh( (a + b * x * x) * x ) + Type(.5);
}
VECTORIZE1(pnorm_approx);

/** \brief  Approximate inverse normal cumulative distribution function, similar to R's qnorm (one-argument case only).
* \details
To be replaced by more accurate version based on Rmath library.
*/
template<class Type>
Type qnorm_approx(Type x){
  Type a=993./880.;
  Type b=89./880.;
  Type y = .5*(log(x)-log(1-x));
  Type p = a/b;
  Type q = -y/b;
  Type Delta0 = -3*p;
  Type Delta1 = 27*q;
  Type C = pow( .5 * Delta1 + .5 * sqrt( pow(Delta1,2) - 4 * pow(Delta0,3) ), Type(1)/Type(3) );
  return -(C + Delta0 / C) * sqrt(Type(2)) / Type(3);
}
VECTORIZE1(qnorm_approx);

/**	\brief Polynomial evaluation.

	Evaluates the given polynomial of degree N at x.
	*/
template <class Type>
Type polevl(Type x, const vector<Type> &coef, int N)
{
	Type res = coef[0];

	for(int i=1;i<=N;i++)
		res = res * x + coef[i];
	
	return(res);
}

#define LOG2E 1.4426950408889634073599 // 1/log(2)

/**	\brief Inverse hyperbolic sine function.

	Based on __Cephes Math Library__ Release 2.8:  June, 2000.
	Copyright 1984, 1995, 2000 by Stephen L. Moshier ; according to http://www.netlib.org/cephes/readme.
	*/
template <class Type>
Type asinh(Type x)
{
	vector<Type> P(5);
	P[0] = Type(-4.33231683752342103572E-3);
	P[1] = Type(-5.91750212056387121207E-1);
	P[2] = Type(-4.37390226194356683570E0);
	P[3] = Type(-9.09030533308377316566E0);
	P[4] = Type(-5.56682227230859640450E0);

	vector<Type> Q(5);
	Q[0] = Type(1);
	Q[1] = Type(1.28757002067426453537E1);
	Q[2] = Type(4.86042483805291788324E1);
	Q[3] = Type(6.95722521337257608734E1);
	Q[4] = Type(3.34009336338516356383E1);

	Type a,z,y,res0,res;
	Type sign;

	if(isnan(x)) return x;
	
	sign = CppAD::CondExpLt(x,Type(0),Type(-1),Type(1));
	y = sign*x;
	
	res0 = CppAD::CondExpEq(y,Type(INFINITY),x,sign*(log(y)+LOG2E));
	res = CppAD::CondExpGt(x,Type(10e8),res0,res);
	
	z = y*y;
	
	a = (polevl(z,P,4)/polevl(z,Q,4))*z;
	a = sign*(a*x+x);
	a = CppAD::CondExpLt(y,Type(0.5),a,sqrt(z+1));
	
	res = CppAD::CondExpLe(x,Type(10e8),sign*log(x+a),res);
	res = CppAD::CondExpLt(y,Type(0.5),a,res);
	res = CppAD::CondExpEq(x,Type(0),x,res);
	
	return sign*res;
}

