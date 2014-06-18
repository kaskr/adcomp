#include <cmath>

/**	\file
	\brief Probability distribution functions.
	*/
	
/**	\brief Cumulative distribution function of the exponential distribution.
	\ingroup R_style_distribution
	\param rate Rate parameter. Must be strictly positive.
	\param give_log 1 if one wants the log-probability, 0 otherwise.
	*/
template<class Type> 
Type pexp(Type x, Type rate, int give_log)
{
	if(!give_log)
	{
		if(x>=0) return 1-exp(-rate*x);
		else return 0;
	} 
	else
	{
		if(x>=0) return log(1-exp(-rate*x));
		else return -INFINITY;
	}
}

// Vectorize pexp
VECTORIZE3_Tti(pexp);
VECTORIZE3_tTi(pexp);
VECTORIZE3_TTi(pexp);

/** 	\brief Cumulative distribution function of the Weibull distribution.
	\ingroup R_style_distribution
	\param shape Shape parameter. Must be strictly positive.
	\param scale Scale parameter. Must be strictly positive.
	\param give_log 1 if one wants the log-probability, 0 otherwise.
	*/
template<class Type> 
Type pweibull(Type x, Type shape, Type scale, int give_log)
{
	if(!give_log)
	{
		if(x>=0) return 1-exp(-pow(x/scale,shape));
		else return 0;
	} 
	else
	{
		if(x>=0) return log(1-exp(-pow(x/scale,shape)));
		else return -INFINITY;
	}
}

// Vectorize pweibull
VECTORIZE4_Ttti(pweibull);
VECTORIZE4_tTti(pweibull);
VECTORIZE4_ttTi(pweibull);
VECTORIZE4_TTti(pweibull);
VECTORIZE4_TtTi(pweibull);
VECTORIZE4_tTTi(pweibull);
VECTORIZE4_TTTi(pweibull);

int fac(int n) // Factorial function.
{
	return (n == 1 || n == 0) ? 1 : fac(n - 1) * n;
}

/**	\brief Probability mass function of the binomial distribution.
	\ingroup R_style_distribution
	\param k Number of successes.
	\param size Number of trials.
	\param prob Probability of success.
	\param give_log 1 if one wants the log-probability, 0 otherwise.
	*/
template<class Type> 
Type dbinom(int k, int size, Type prob, int give_log)
{
	if(!give_log) return Type(fac(size))/Type(fac(k)*fac(size-k))*pow(prob,k)*pow(1-prob,size-k);
	else return log(fac(size))-log(fac(k))-log(fac(size-k))+k*log(prob)+(size-k)*log(1-prob);
}

// Vectorize dbinom
VECTORIZE4_Iiti(dbinom);
VECTORIZE4_iIti(dbinom);
VECTORIZE4_iiTi(dbinom);
VECTORIZE4_IIti(dbinom);
VECTORIZE4_IiTi(dbinom);
VECTORIZE4_iITi(dbinom);
VECTORIZE4_IITi(dbinom);

/**	\brief Probability density function of the exponential distribution.
	\ingroup R_style_distribution
	\param rate Rate parameter. Must be strictly positive.
	\param give_log 1 if one wants the log-probability, 0 otherwise.
	*/
template<class Type> 
Type dexp(Type x, Type rate, int give_log)
{
	if(!give_log)
	{
		if(x>=0) return rate*exp(-rate*x);
		else return 0;
	} 
	else
	{
		if(x>=0) return log(rate) - rate*x;
		else return -INFINITY;
	}
}

// Vectorize dexp
VECTORIZE3_Tti(dexp);
VECTORIZE3_tTi(dexp);
VECTORIZE3_TTi(dexp);

/** 	\brief Probability density function of the Weibull distribution.
	\ingroup R_style_distribution
	\param shape Shape parameter. Must be strictly positive.
	\param scale Scale parameter. Must be strictly positive.
	\param give_log 1 if one wants the log-probability, 0 otherwise.
	*/
template<class Type> 
Type dweibull(Type x, Type shape, Type scale, int give_log)
{
	if(!give_log)
	{
		if(x>=0) return shape/scale * pow(x/scale,shape-1) * exp(-pow(x/scale,shape));
		else return 0;
	} 
	else
	{
		if(x>=0) return log(shape) - log(scale) + (shape-1)*(log(x)-log(scale)) - pow(x/scale,shape);
		else return -INFINITY;
	}
}

// Vectorize dweibull
VECTORIZE4_Ttti(dweibull);
VECTORIZE4_tTti(dweibull);
VECTORIZE4_ttTi(dweibull);
VECTORIZE4_TTti(dweibull);
VECTORIZE4_TtTi(dweibull);
VECTORIZE4_tTTi(dweibull);
VECTORIZE4_TTTi(dweibull);

/**	\brief Inverse cumulative distribution function of the Weibull distribution.
	\ingroup R_style_distribution
	\param p Probability ; must be between 0 and 1.
	\param shape Shape parameter. Must be strictly positive.
	\param scale Scale parameter. Must be strictly positive.
	\param give_log 1 if one wants the log-probability, 0 otherwise.
	*/
template<class Type> 
Type qweibull(Type p, Type shape, Type scale, int give_log)
{
	if(p<0 || p>1) return 0;
	if(!give_log) return scale * pow( (-log(1-p)) , 1/shape );
	else return log(scale) + 1/shape * log(-log(1-p));
}

// Vectorize qweibull
VECTORIZE4_Ttti(qweibull);
VECTORIZE4_tTti(qweibull);
VECTORIZE4_ttTi(qweibull);
VECTORIZE4_TTti(qweibull);
VECTORIZE4_TtTi(qweibull);
VECTORIZE4_tTTi(qweibull);
VECTORIZE4_TTTi(qweibull);
