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
Type polevl(Type x, Type *coef, int N)
{
	Type *p = coef;
	Type ans = *p;
 
	p = p+1;

	for(int i=1;i<=N;i++)
	{
		ans = ans * x + *p;
		p = p+1;
	}
	
	return( ans );
}

/**	\brief Polynomial evaluation.

	Evaluates the given polynomial of degree N at x, assuming coefficient of N is 1.0.
	*/
template <class Type> 
Type p1evl(Type x, Type *coef, int N)
{
	Type *p = coef;
	Type ans = x + *p;
	
	p = p+1;

	for(int i=1;i<=N-1;i++)
	{
		ans = ans * x + *p;
		p = p+1;
	}

	return( ans );
}

#define LOG2E 1.4426950408889634073599 // 1/log(2)

/**	\brief Inverse hyperbolic sine function.
	*/
template <class Type>
Type asinh(Type x)
{
	Type P[] = {
		-4.33231683752342103572E-3,
		-5.91750212056387121207E-1,
		-4.37390226194356683570E0,
		-9.09030533308377316566E0,
		-5.56682227230859640450E0
	};
	Type Q[] = {
		1.28757002067426453537E1,
		4.86042483805291788324E1,
		6.95722521337257608734E1,
		3.34009336338516356383E1
	};

	Type a,z,y;
	int sign;

	if(isnan(x)) return x;
	if(x==Type(0)) return x;
	
	if(x<Type(0))
	{
		sign = -1;
		y = -x;
	}
	else
	{
		sign = 1;
		y = x;
	}
	
	if(x>Type(10e8))
	{
		if(y==INFINITY) return x;
		return sign*(log(x)+LOG2E);
	}
	
	z = y*y;
	
	if(y<Type(0.5))
	{
		a = (polevl(z,P,4)/p1evl(z,Q,4))*z;
		a = a*x+x;
		if(sign<0) a = -a;
		return a;
	}
	
	a = sqrt(z+1.);
	return sign*log(x+a);
}

