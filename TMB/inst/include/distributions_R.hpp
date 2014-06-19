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
		return CppAD::CondExpGe(x,Type(0),1-exp(-rate*x),Type(0));
	else
		return CppAD::CondExpGe(x,Type(0),log(1-exp(-rate*x)),Type(-INFINITY));
}

// Vectorize pexp
VECTORIZE3_tti(pexp);

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
		return CppAD::CondExpGe(x,Type(0),1-exp(-pow(x/scale,shape)),Type(0));
	else
		return CppAD::CondExpGe(x,Type(0),log(1-exp(-pow(x/scale,shape))),Type(-INFINITY));
}

// Vectorize pweibull
VECTORIZE4_ttti(pweibull);

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
	Type logres = lgamma(size+1)-lgamma(k+1)-lgamma(size-k+1)+k*log(prob)+(size-k)*log(1-prob);
	if(!give_log) return exp(logres);
	else return logres;
}

// Vectorize dbinom
VECTORIZE4_iiti(dbinom);

/**	\brief Probability density function of the exponential distribution.
	\ingroup R_style_distribution
	\param rate Rate parameter. Must be strictly positive.
	\param give_log 1 if one wants the log-probability, 0 otherwise.
	*/
template<class Type> 
Type dexp(Type x, Type rate, int give_log)
{
	if(!give_log)
		return CppAD::CondExpGe(x,Type(0),rate*exp(-rate*x),Type(0));
	else
		return CppAD::CondExpGe(x,Type(0),log(rate)-rate*x,Type(-INFINITY));
}

// Vectorize dexp
VECTORIZE3_tti(dexp);

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
		return CppAD::CondExpGe(x,Type(0),shape/scale * pow(x/scale,shape-1) * exp(-pow(x/scale,shape)),Type(0));
	else
		return CppAD::CondExpGe(x,Type(0),log(shape) - log(scale) + (shape-1)*(log(x)-log(scale)) - pow(x/scale,shape),Type(-INFINITY));
}

// Vectorize dweibull
VECTORIZE4_ttti(dweibull);

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
	Type res;	
	
	if(!give_log) res = scale * pow( (-log(1-p)) , 1/shape );
	else res = log(scale) + 1/shape * log(-log(1-p));
	
	res = CppAD::CondExpLt(p,Type(0),Type(0),res);
	res = CppAD::CondExpGt(p,Type(1),Type(0),res);
	return res;
}

// Vectorize qweibull
VECTORIZE4_ttti(qweibull);


