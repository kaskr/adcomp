// Copyright (C) 2013-2015 Kasper Kristensen
// License: GPL-2

/** \file
    \brief Univariate normal density
    \ingroup R_style_distribution
*/

/** \brief Probability density function of the normal distribution
    \ingroup R_style_distribution
    \param x vector of observations
    \param mean mean of the normal distribution
    \param sd standard deviation of the normal distribution. must be strictly positive.
    \param give_log true if one wants the log-probability, false otherwise.
*/
template<class Type>
Type dnorm(Type x, Type mean, Type sd, int give_log=0)
{
  Type resid = (x - mean) / sd;
  Type logans = Type(-log(sqrt(2*M_PI))) - log(sd) - Type(.5) * resid * resid;
  if(give_log) return logans; else return exp(logans);
}
VECTORIZE4_ttti(dnorm)
