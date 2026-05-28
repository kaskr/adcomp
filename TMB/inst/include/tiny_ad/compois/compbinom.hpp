namespace compbinom_utils {

/** \brief Conway-Maxwell-Binomial. Calculate log-normalizing constant.
 *
 *  logZ(Psi, nu, n) = log sum_{k=0}^n C(n,k)^nu * exp(k*Psi)
 *
 *  Direct summation over n+1 terms in log space. No truncation needed
 *  since the support is finite. log C(n,k) is built incrementally via
 *  the standard binomial-coefficient recursion for numerical stability.
 *
 *  \param Psi  log(p/(1-p)), the natural (logit) parameter
 *  \param nu   dispersion parameter (>0)
 *  \param n    number of trials (treated as constant, never differentiated)
 *  \return     log Z
 */
template<class Type>
Type calc_logZ(Type Psi, Type nu, int n) {
  using atomic::tiny_ad::isfinite;
  bool ok = (n >= 0 && isfinite(Psi) && isfinite(nu));
  if (!ok) return NAN;
  using atomic::robust_utils::logspace_add;
  Type logZ = Type(-INFINITY);
  Type log_choose_k = Type(0);  // log C(n, 0) = 0
  for (int k = 0; k <= n; ++k) {
    if (k > 0) {
      // log C(n, k) = log C(n, k-1) + log(n-k+1) - log(k)
      log_choose_k += log((double)(n - k + 1)) - log((double)k);
    }
    Type log_term = nu * log_choose_k + (double)k * Psi;
    logZ = logspace_add(logZ, log_term);
  }
  return logZ;
}

/** \brief Conway-Maxwell-Binomial. Calculate mean from natural parameter. */
template<class Type>
Type calc_mean(Type Psi, Type nu, int n) {
  typedef atomic::tiny_ad::variable<1, 1, Type> ADType;
  ADType Psi_(Psi, 0);
  ADType ans = calc_logZ<ADType>(Psi_, nu, n);
  return ans.getDeriv()[0];
}

/** \brief Conway-Maxwell-Binomial. Calculate logit(p) from log(mean).
 *
 *  Inverts E[Y | n, Psi, nu] = exp(log_mean) for Psi using safeguarded
 *  Newton iteration with bisection backstop. Pure Newton oscillates at
 *  strong overdispersion (nu << 1); the bracket-based variant falls back
 *  to bisection when the Newton step lands outside [Psi_lo, Psi_hi].
 */
template<class Type>
Type calc_logitp(Type log_mean, Type nu, int n) {
  using atomic::tiny_ad::isfinite;
  bool ok = (n >= 0 && isfinite(log_mean) && isfinite(nu));
  if (!ok) return NAN;
  int iter_max = 200;
  double reltol = 1e-12;
  double abstol = 1e-14;
  typedef atomic::tiny_ad::variable<1, 1, Type> ADType;
  Type mu = exp(log_mean);
  Type Psi_lo = Type(-30.0);
  Type Psi_hi = Type(+30.0);
  ADType x(log(mu / (Type(n) - mu)), 0);
  int i;
  for (i = 0; i < iter_max; i++) {
    x.deriv[0] = Type(1.0);
    ADType y = calc_mean<ADType>(x, nu, n);
    Type residual = y.value - mu;
    if (residual > Type(0)) {
      Psi_hi = x.value;
    } else {
      Psi_lo = x.value;
    }
    if (fabs(residual) <= reltol * fabs(mu)) break;
    if (fabs(residual) <= abstol) break;
    Type step;
    if (y.deriv[0] > Type(1e-300)) {
      step = -residual / y.deriv[0];
    } else {
      step = (Psi_lo + Psi_hi) / Type(2) - x.value;
    }
    Type x_next = x.value + step;
    if (x_next > Psi_lo && x_next < Psi_hi) {
      x.value = x_next;
    } else {
      x.value = (Psi_lo + Psi_hi) / Type(2);
    }
    if ((Psi_hi - Psi_lo) < Type(abstol)) break;
  }
  if (i == iter_max) {
    Rf_warning("compbinom_utils::calc_logitp: Maximum number of iterations exceeded");
  }
  return x.value;
}

} // namespace compbinom_utils
