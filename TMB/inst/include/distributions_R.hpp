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
	\param give_log 1 if one wants the log-probability, 0 otherwise.
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
	\param give_log 1 if one wants the log-probability, 0 otherwise.
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
/**@}*/


/**	@name Weibull distribution.
	Functions relative to the Weibull distribution.
	*/
/**@{*/
/** 	\brief Cumulative distribution function of the Weibull distribution.
	\ingroup R_style_distribution
	\param shape Shape parameter. Must be strictly positive.
	\param scale Scale parameter. Must be strictly positive.
	\param give_log 1 if one wants the log-probability, 0 otherwise.
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
	\param give_log 1 if one wants the log-probability, 0 otherwise.
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
	\param give_log 1 if one wants the log-probability, 0 otherwise.
	*/
template<class Type> 
Type qweibull(Type p, Type shape, Type scale, int give_log=0)
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
/**@}*/

/**	\brief Probability mass function of the binomial distribution.
	\ingroup R_style_distribution
	\param k Number of successes.
	\param size Number of trials.
	\param prob Probability of success.
	\param give_log 1 if one wants the log-probability, 0 otherwise.
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
	\param give_log 1 if one wants the log-probability, 0 otherwise.
	*/
template <class Type>
Type dbeta(Type x, Type shape1, Type shape2, int give_log)
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
	\param give_log 1 if one wants the log-probability, 0 otherwise.
	*/
template <class Type>
Type df(Type x, Type df1, Type df2, int give_log)
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
	\param give_log 1 if one wants the log-probability, 0 otherwise.
	*/
template <class Type>
Type dlogis(Type x, Type location, Type scale, int give_log)
{
	Type logres = -(x-location)/scale - log(scale) - 2*log(1+exp(-(x-location)/scale));
	if(!give_log) return exp(logres);
	else return logres;
}

// Vectorize dlogis
VECTORIZE4_ttti(dlogis);

/**	@name Normal distribution.
	Functions relative to the normal distribution.
	*/
/**@{*/
/**	\brief Cumulative distribution function of the standard normal distribution.
	\ingroup R_style_distribution
	\param give_log 1 if one wants the log-probability, 0 otherwise.
	*/
template<class Type>
Type pnorm_standard(Type x, int give_log)
{
	return Type(0.5)+0.5*erf(x/sqrt(2));
}

/**	\brief Cumulative distribution function of the normal distribution.
	\ingroup R_style_distribution
	\param mean Mean of the normal distribution.
	\param sd Standard deviation of the normal distribution.
	\param give_log 1 if one wants the log-probability, 0 otherwise.
	*/
template<class Type>
Type pnorm(Type x, Type mean, Type sd, int give_log)
{
	Type p = (x-mean)/sd;

	Type tmp = CppAD::CondExpLt(x,mean,Type(0),Type(1));
	return CppAD::CondExpLe(sd,Type(0),tmp,pnorm_standard(p,give_log));
}

// Vectorize pnorm
VECTORIZE4_ttti(pnorm);
/**@}*/

/**	\brief Probability density function of the skew-normal distribution.
	\ingroup R_style_distribution
	\param alpha Slant parameter.
	\param give_log 1 if one wants the log-probability, 0 otherwise.
	*/
template <class Type>
Type dsn(Type x, Type alpha, int give_log=0)
{
	if(!give_log) return 2 * dnorm<Type>(x,Type(0),Type(1),0) * pnorm<Type>(alpha*x,Type(0),Type(1),0);
	else return log(2) + log(dnorm<Type>(x,Type(0),Type(1),0)) + log(pnorm<Type>(alpha*x,Type(0),Type(1),0));
}

// Vectorize dsn
VECTORIZE3_tti(dsn);

/**	@name Tweedie distribution.
	Functions relative to the Tweedie distribution.
	*/
/**@{*/
/**	\brief Deviance function for the Tweedie distribution.
	\ingroup R_style_distribution
	\param mu Mean.
	\param power Value of p such that the variance is \f$ var[Y] = \phi*\mu^p \f$.
	*/
template <class Type>
Type tweedie_dev(Type y, Type mu, Type power)
{
	Type p = power;
	Type res;
	Type dev0, dev1, dev2, dev3;
	
	dev0 = (y-mu)*(y-mu)/2;
	dev1 = CppAD::CondExpEq(y,Type(0),mu,y*log(y/mu)-(y-mu));
	dev2 = log(mu/y)+y/mu-1;
	dev3 = pow(y,2-p)/((1-p)*(2-p)) - (y*pow(mu,1-p))/(1-p) + pow(mu,2-p)/(2-p);
	res = dev3;

	res = CppAD::CondExpEq(p,Type(0),dev0,res);
	res = CppAD::CondExpEq(p,Type(1),dev1,res);
	res = CppAD::CondExpEq(p,Type(2),dev2,res);
	return 2*res;
}

/**	\brief Saddlepoint density for the Tweedie distribution.
	\ingroup R_style_distribution
 	\param xi Value of p such that the variance is \f$ var[Y] = \phi*\mu^p \f$.
	\param mu Mean.
	\param phi Dispersion.
	\param eps Offset in computing the variance function. Default value : 1/6.
	\param give_log 1 if one wants the log-probability, 0 otherwise (default).
	*/	
template <class Type>
Type dtweedie_saddle(Type y, Type xi, Type mu, Type phi, Type eps=Type(1)/6, int give_log=0)
{
	Type p = xi;
	Type y_eps = y;
	Type dev;
	Type density0, density;
	
		
	y_eps = CppAD::CondExpLt(p,Type(2),y+eps,y_eps);
	
	dev = tweedie_dev(y,mu,p);
	density = 1/sqrt(2*M_PI*phi*pow(y_eps,p))*exp(-dev/(2*phi));

	density0 = CppAD::CondExpLt(p,Type(2), exp(-pow(mu,2-p)/(phi*(2-p))) ,Type(0));
	density0 = CppAD::CondExpLt(p,Type(1),Type(0),density0);
	
	density = CppAD::CondExpEq(y,Type(0),density0,density);
	if(!give_log) return density;
	else return log(density);
}

// Vectorize dtweedie_saddle
VECTORIZE6_Ttttti(dtweedie_saddle);


/**@}*/

