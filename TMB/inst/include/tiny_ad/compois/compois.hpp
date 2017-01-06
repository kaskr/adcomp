namespace compois_utils {
//#define TRACE(x) std::cout << #x << "=" << x << "\n";

/** \brief Conway-Maxwell-Poisson. Calculate log-normalizing
    constant. */
template<class Type>
Type calc_logZ(Type loglambda, Type nu){
  using atomic::tiny_ad::isfinite;
  bool ok = (0 < nu && isfinite(loglambda) && isfinite(nu) );
  if (!ok) return NAN;
  using atomic::robust_utils::logspace_add;
  using atomic::robust_utils::logspace_sub;
  int iter_max = 1e4;
  Type logZ = 0.;
  Type logmu = loglambda / nu;
  Type mu = exp(logmu); // mode: lambda^(1/nu)
  if (false) {
    // Asymptotic expansion for large mu (dropped because inaccurate)
    logZ =
      nu * mu - 0.5 * log(nu) -
      ((nu - 1) / 2) * log(mu) -
      ((nu - 1) / 2) * log(2 * M_PI);
  } else if ( (mu > 100) && (mu * nu > 200) && (nu < 2 * mu) ) {
    using atomic::tiny_ad::lgamma;
    // Laplace approximation when nu=1 (Poisson case)
    Type jhat = mu - .5;
    Type H1 = lgamma<2>(jhat + 1.);
    Type fhat1 = jhat * logmu - lgamma(jhat + 1.);
    Type logL1 = log(sqrt(2. * M_PI)) - .5 * log(H1) + fhat1;
    // Laplace approximation error nu=1
    Type err1 = logL1 - mu;
    // Laplace approximation general nu
    Type H = nu * H1;
    Type fhat = nu * fhat1;
    Type logL = log(sqrt(2. * M_PI)) - .5 * log(H) + fhat;
    // Apply correction
    logL -= err1 / nu;
    logZ = logL;
  }
  else {
    // Series summation
    int index;             // Current index
    Type logT;             // Current log term
    Type dlogT;            // logT(index) - logT(index-1)
    double reltol = 1e-12; // Break if Term_current / Sum_current < reltol
    int i;
    // Initialize largest term and sum
    int index_mode = floor(mu);
    Type logT_mode = index_mode * loglambda - nu * lgamma(index_mode + 1.);
    logZ = logT_mode;
    // Left tail. FIXME: Use half gaussian approximation for mu >~ maxit.
    logT = logT_mode; // Initialize
    for(i = 1; i<iter_max; i++) {
      index = index_mode - i;
      if (index < 0) break;
      dlogT = loglambda - nu * log( (double) index + 1);
      logT -= dlogT;
      logZ = logspace_add(logZ, logT);
      if ( logT - logZ < log(reltol) ) break;
    }
    // Right tail
    logT = logT_mode; // Initialize
    for(i = 1; i<iter_max; i++) {
      index = index_mode + i;
      dlogT = loglambda - nu * log( (double) index);
      logT += dlogT;
      logZ = logspace_add(logZ, logT);
      if ( logT - logZ < log(reltol) ) break;
    }
    // Tail upper bound via geometric series
    // T_j = T_i * exp( dlogT_i * j ) , j>=i
    //  sum( rho^j , j=i,...) = rho^i * 1/(1-rho)
    // T_i * rho^i * 1/(1-rho)
    // logT_tail = log(T_i) + i * log(rho) - log(1-rho)
    Type logT_tail = logT + index * dlogT - logspace_sub(Type(0), dlogT);
    logZ = logspace_add(logZ, logT_tail);
  }
  return logZ;
}

/** \brief Conway-Maxwell-Poisson. Calculate mean from log(lambda). */
template<class Type>
Type calc_mean(Type loglambda, Type nu){
  typedef atomic::tiny_ad::variable<1, 1, Type> ADType;
  ADType loglambda_ (loglambda, 0);
  ADType ans = calc_logZ<ADType>(loglambda_, nu);
  return ans.getDeriv()[0];
}

/** \brief Conway-Maxwell-Poisson. Calculate log(lambda) from
    log(mean). */
template<class Type>
Type calc_loglambda(Type logmean, Type nu) {
  using atomic::tiny_ad::isfinite;
  bool ok = (0 < nu && isfinite(logmean) && isfinite(nu) );
  if (!ok) return NAN;
  int iter_max = 100; double reltol = 1e-9, abstol = 1e-12;
  typedef atomic::tiny_ad::variable<1, 1, Type> ADType;
  ADType x = nu * logmean; // Initial guess
  ADType step = 0;
  ADType f_previous = INFINITY;
  int i=0;
  for ( ; i < iter_max ; i++ ) {
    x.deriv[0] = 1.; // Seed
    ADType y = calc_mean<ADType>(x, nu);
    if( ! isfinite(y) ) {
      if (i==0) return NAN; // Undefined initial value
      // Step halfway back
      step = step / 2;
      x -= step;
      continue;
    }
    ADType f = ( y > 0 ?
                 log(y) - ADType( logmean ) :
                 y - ADType( exp(logmean) ) );
    if( fabs(f) > fabs(f_previous) ) {
      // Step halfway back
      step = step / 2;
      x -= step;
      continue;
    }
    step = ( f.deriv[0] != 0 ?
             -f.value / f.deriv[0] : 0 ); // Newton step
    x += step; f_previous = f;
    if(fabs(step) <= reltol * fabs(x))
      break;
    if(fabs(step) <= abstol)
      break;
  }
  if (i == iter_max)
    warning("calc_loglambda: Maximum number of iterations exceeded");
  return x.value;
}

} // compois_utils
