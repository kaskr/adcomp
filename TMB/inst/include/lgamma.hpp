template <class Type>
inline Type lgamma(const Type &y)
{
  /* coefficients for gamma=7, kmax=8  Lanczos method */
  static const Type
    LogRootTwoPi_ = 0.9189385332046727418,
    lanczos_7_c[9] = {
      0.99999999999980993227684700473478,
      676.520368121885098567009190444019,
      -1259.13921672240287047156078755283,
      771.3234287776530788486528258894,
      -176.61502916214059906584551354,
      12.507343278686904814458936853,
      -0.13857109526572011689554707,
      9.984369578019570859563e-6,
      1.50563273514931155834e-7
    };
  Type x=y;
  int k;
  Type Ag;
  Type term1, term2;
  x -= Type(1.0); /* Lanczos writes z! instead of Gamma(z) */
  Ag = lanczos_7_c[0];
  for(k=1; k<=8; k++) { Ag += lanczos_7_c[k]/(x+k); }
  /* (x+0.5)*log(x+7.5) - (x+7.5) + LogRootTwoPi_ + log(Ag(x)) */
  term1 = (x+Type(0.5))*log((x+Type(7.5))/Type(M_E));
  term2 = LogRootTwoPi_ + log(Ag);
  return term1 + (term2 - Type(7.0));
}

/* Negative binomial probability function.
 * Considering the return value we need to make sure that:
 * (1) n>0   (2) 0<p   (3) p<1
 * This is obtained by adding small constants appropriate places.
 */
template<class Type>
inline Type dnbinom(const Type &x, const Type &mu0, const Type &var0,
		    int give_log=0)
{
  Type mu=mu0+Type(1e-15);
  Type var=var0+Type(2e-15);
  Type p=mu/var;
  Type n=mu*p/(1-p)+Type(1e-15);
  Type logres = lgamma(x+n)-lgamma(n)-lgamma(x+Type(1))+
    n*log(p)+x*log(Type(1)-p);
  if (give_log) return logres; else return exp(logres);
}

/* Poisson probability function. */
template<class Type>
inline Type dpois(const Type &x, const Type &lambda, int give_log=0)
{
  Type logres = -lambda + x*log(lambda) - lgamma(x+Type(1));
  if (give_log) return logres; else return exp(logres);
}

/* Density of X where X~gamma distributed */
template<class Type>
Type dgamma(Type y, Type shape, Type scale, int give_log=0)
{
  Type logres=-lgamma(shape)+(shape-Type(1.0))*log(y)-y;
  if(give_log)return logres; else return exp(logres);
}

/* Density of log(X) where X~gamma distributed */
template<class Type>
inline Type dlgamma(Type y, Type shape, Type scale, int give_log=0)
{
  Type logres=-lgamma(shape)-shape*log(scale)-exp(y)/scale+shape*y;
  if(give_log)return logres; else return exp(logres);
}
