/**	\file
	\brief Probability distribution functions.
	*/

/** 	@name Exponential distribution.
	Functions relative to the exponential distribution.
	*/
/**@{*/
/**	\brief Cumulative distribution function of the exponential distribution.
	\ingroup R_style_distribution
	\param rate Rate parameter. Must be strictly positive.
	\param give_log true if one wants the log-probability, false otherwise.
	*/
template<class Type> 
Type pexp(Type x, Type rate, int give_log=0)
{
	if(!give_log)
		return CppAD::CondExpGe(x,Type(0),1-exp(-rate*x),Type(0));
	else
		return CppAD::CondExpGe(x,Type(0),log(1-exp(-rate*x)),Type(-INFINITY));
}

// Vectorize pexp
VECTORIZE3_tti(pexp);

/**	\brief Probability density function of the exponential distribution.
	\ingroup R_style_distribution
	\param rate Rate parameter. Must be strictly positive.
	\param give_log true if one wants the log-probability, false otherwise.
	*/
template<class Type> 
Type dexp(Type x, Type rate, int give_log=0)
{
	if(!give_log)
		return CppAD::CondExpGe(x,Type(0),rate*exp(-rate*x),Type(0));
	else
		return CppAD::CondExpGe(x,Type(0),log(rate)-rate*x,Type(-INFINITY));
}

// Vectorize dexp
VECTORIZE3_tti(dexp);

/**	\brief Inverse cumulative distribution function of the exponential distribution.
	\ingroup R_style_distribution
	\param rate Rate parameter. Must be strictly positive.
	\param log_p true if p is log-probability, false otherwise.
	*/
template <class Type>
Type qexp(Type p, Type rate, int log_p=0)
{
	if(!log_p) return -log(1-p)/rate;
	else return -log(1-exp(p))/rate;
}

// Vectorize qexp.
VECTORIZE3_tti(qexp);
/**@}*/


/**	@name Weibull distribution.
	Functions relative to the Weibull distribution.
	*/
/**@{*/
/** 	\brief Cumulative distribution function of the Weibull distribution.
	\ingroup R_style_distribution
	\param shape Shape parameter. Must be strictly positive.
	\param scale Scale parameter. Must be strictly positive.
	\param give_log true if one wants the log-probability, false otherwise.
	*/
template<class Type> 
Type pweibull(Type x, Type shape, Type scale, int give_log=0)
{
	if(!give_log)
		return CppAD::CondExpGe(x,Type(0),1-exp(-pow(x/scale,shape)),Type(0));
	else
		return CppAD::CondExpGe(x,Type(0),log(1-exp(-pow(x/scale,shape))),Type(-INFINITY));
}

// Vectorize pweibull
VECTORIZE4_ttti(pweibull);

/** 	\brief Probability density function of the Weibull distribution.
	\ingroup R_style_distribution
	\param shape Shape parameter. Must be strictly positive.
	\param scale Scale parameter. Must be strictly positive.
	\param give_log true if one wants the log-probability, false otherwise.
	*/
template<class Type> 
Type dweibull(Type x, Type shape, Type scale, int give_log=0)
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
	\param log_p true if p is log-probability, false otherwise.
	*/
template<class Type> 
Type qweibull(Type p, Type shape, Type scale, int log_p=0)
{
	Type res;
	
	if(!log_p) res = scale * pow( (-log(1-p)) , 1/shape );
	else res = scale * pow( (-log(1-exp(p))) , 1/shape );
	res = CppAD::CondExpLt(p,Type(0),Type(0),res);
	res = CppAD::CondExpGt(p,Type(1),Type(0),res);
	return res;
}

// Vectorize qweibull
VECTORIZE4_ttti(qweibull);
/**@}*/

/**	\brief Probability mass function of the binomial distribution.
	\ingroup R_style_distribution
	\param k Number of successes.
	\param size Number of trials.
	\param prob Probability of success.
	\param give_log true if one wants the log-probability, false otherwise.
	*/
template<class Type> 
Type dbinom(Type k, Type size, Type prob, int give_log=0)
{
	Type logres = lgamma(size+1)-lgamma(k+1)-lgamma(size-k+1)+k*log(prob)+(size-k)*log(1-prob);
	if(!give_log) return exp(logres);
	else return logres;
}

// Vectorize dbinom
VECTORIZE4_ttti(dbinom);

/**	\brief Probability density function of the beta distribution.
	\ingroup R_style_distribution
	\param shape1 First shape parameter. Must be strictly positive.
	\param shape2 Second shape parameter. Must be strictly positive.
	\param give_log true if one wants the log-probability, false otherwise.
	*/
template <class Type>
Type dbeta(Type x, Type shape1, Type shape2, int give_log=0)
{
	Type res = exp(lgamma(shape1+shape2) - lgamma(shape1) - lgamma(shape2)) * pow(x,shape1-1) * pow(1-x,shape2-1);
	if(!give_log) 
		return res;
	else 
		return CppAD::CondExpEq(x,Type(0),log(res),lgamma(shape1+shape2) - lgamma(shape1) - lgamma(shape2) + (shape1-1)*log(x) + (shape2-1)*log(1-x));
}

// Vectorize dbeta
VECTORIZE4_ttti(dbeta);

/**	\brief Probability density function of the Fisher distribution.
	\ingroup R_style_distribution
	\param df1 Degrees of freedom 1.
	\param df2 Degrees of freedom 2.
	\param give_log true if one wants the log-probability, false otherwise.
	*/
template <class Type>
Type df(Type x, Type df1, Type df2, int give_log=0)
{
	Type logres = lgamma((df1+df2)/2.) - lgamma(df1/2.) - lgamma(df2/2.) + df1/2.*log(Type(df1)/df2) + (df1/2.-1)*log(x) - (df1+df2)/2.*log(1+Type(df1)/df2*x);
	if(!give_log) return exp(logres);
	else return logres;
}

//Vectorize df
VECTORIZE4_ttti(df);

/**	\brief Probability density function of the logistic distribution.
	\ingroup R_style_distribution
	\param location Location parameter.
	\param scale Scale parameter. Must be strictly positive.
	\param give_log true if one wants the log-probability, false otherwise.
	*/
template <class Type>
Type dlogis(Type x, Type location, Type scale, int give_log=0)
{
	Type logres = -(x-location)/scale - log(scale) - 2*log(1+exp(-(x-location)/scale));
	if(!give_log) return exp(logres);
	else return logres;
}

// Vectorize dlogis
VECTORIZE4_ttti(dlogis);

/**	\brief Probability density function of the skew-normal distribution.
	\ingroup R_style_distribution
	\param alpha Slant parameter.
	\param give_log true if one wants the log-probability, false otherwise.
	
	Notation adopted from R package "sn".

	The skew normal distribution generalises the normal distribution to allow for non-zero skewness.
	*/
template <class Type>
Type dsn(Type x, Type alpha, int give_log=0)
{
	// TODO : change pnorm_approx to pnorm when pnorm is written	
	
	if(!give_log) return 2 * dnorm(x,Type(0),Type(1),0) * pnorm_approx(alpha*x);
	else return log(2) + log(dnorm(x,Type(0),Type(1),0)) + log(pnorm_approx(alpha*x));
}

// Vectorize dsn
VECTORIZE3_tti(dsn);

/** 	\brief Probability density function of the Student t-distribution.
	\ingroup R_style_distribution
	\param df Degree of freedom.
	\param give_log true if one wants the log-probability, false otherwise.
	*/	
template <class Type>
Type dt(Type x, Type df, int give_log = 0)
{
	Type logres = lgamma((df+1)/2) - Type(1)/2*log(df*M_PI) -lgamma(df/2) - (df+1)/2*log(1+x*x/df);
	if(!give_log) return exp(logres);
	else return logres;
}

// Vectorize dt
VECTORIZE3_tti(dt);

/** 	@name Sinh-asinh distribution.
	Functions relative to the sinh-asinh distribution.
	*/
/**@{*/
/**	\brief Probability density function of the sinh-asinh distribution.
	\ingroup R_style_distribution
	\param mu Location.
	\param sigma Scale.
	\param nu Skewness.
	\param tau Kurtosis.
	\param give_log true if one wants the log-probability, false otherwise.
	
	Notation adopted from R package "gamlss.dist".
	
	Probability density given in (2) in __Jones and Pewsey (2009) Biometrika (2009) 96 (4): 761-780__.
	
	It is not possible to call this function with nu a vector or tau a vector.
	*/
template <class Type>
Type dSHASHo(Type x, Type mu, Type sigma, Type nu, Type tau, int give_log = 0)
{
	Type z = (x-mu)/sigma;
   	Type c = cosh(tau*asinh(z)-nu);
   	Type r = sinh(tau*asinh(z)-nu);
   	Type logres = -log(sigma) + log(tau) -0.5*log(2*M_PI) -0.5*log(1+(z*z)) +log(c) -0.5*(r*r);
   	
   	if(!give_log) return exp(logres);
   	else return logres;
}

// Vectorize dSHASHo
VECTORIZE6_ttttti(dSHASHo);

/**	\brief Cumulative distribution function of the sinh-asinh distribution.
	\ingroup R_style_distribution
	\param mu Location.
	\param sigma Scale.
	\param nu Skewness.
	\param tau Kurtosis.
	\param give_log true if one wants the log-probability, false otherwise.
	
	Notation adopted from R package "gamlss.dist".
	
	It is not possible to call this function with nu a vector or tau a vector.
	*/
template <class Type>
Type pSHASHo(Type q,Type mu,Type sigma,Type nu,Type tau,int give_log=0)
{
	// TODO : Replace pnorm_approx by pnorm when it is written.

  	Type z = (q-mu)/sigma;
  	Type r = sinh(tau * asinh(z) - nu);
  	Type p = pnorm_approx(r);
  	
  	if (!give_log) 
        	return p;
  	else 
  		return log(p);
}

// Vectorize pSHASHo
VECTORIZE6_ttttti(pSHASHo);

/**	\brief Quantile function of the sinh-asinh distribution.
	\ingroup R_style_distribution
	\param mu Location.
	\param sigma Scale.
	\param nu Skewness.
	\param tau Kurtosis.
	\param log_p true if p is log-probability, false otherwise.
	
	Notation adopted from R package "gamlss.dist".
	
	It is not possible to call this function with nu a vector or tau a vector.
	*/
template <class Type>
Type qSHASHo(Type p, Type mu, Type sigma, Type nu, Type tau, int log_p = 0)
{
	// TODO : Replace qnorm_approx by qnorm when it is written.

   	if(!log_p) return mu + sigma*sinh((1/tau)*asinh(qnorm_approx(p))+(nu/tau));
   	else return mu + sigma*sinh((1/tau)*asinh(qnorm_approx(exp(p)))+(nu/tau));
}

// Vectorize qSHASHo
VECTORIZE6_ttttti(qSHASHo);

/**	\brief Transforms a normal variable into a sinh-asinh variable.
	\param mu Location parameter of the result sinh-asinh distribution.
	\param sigma Scale parameter of the result sinh-asinh distribution.
	\param nu Skewness parameter of the result sinh-asinh distribution.
	\param tau Kurtosis parameter of the result sinh-asinh distribution.
	\param log_p true if p is log-probability, false otherwise.
	
	It is not possible to call this function with nu a vector or tau a vector.
	*/
template <class Type>
Type norm2SHASHo(Type x, Type mu, Type sigma, Type nu, Type tau, int log_p = 0)
{
	// TODO : Replace pnorm_approx by pnorm when it is written.

	return qSHASHo(pnorm_approx(x),mu,sigma,nu,tau,log_p);
}

// Vectorize norm2SHASHo
VECTORIZE6_ttttti(norm2SHASHo);
/**@}*/

