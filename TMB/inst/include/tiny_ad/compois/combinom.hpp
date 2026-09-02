namespace combinom_utils {

/** \brief Conway-Maxwell-Binomial. Calculate log-normalizing constant.
 *
 *  logZ(logitp, nu, n) = log sum_{k=0}^n C(n,k)^nu * exp(k*logitp)
 *
 *  Direct summation over n+1 terms in log space. No truncation needed
 *  since the support is finite. log C(n,k) is built incrementally via
 *  the standard binomial-coefficient recursion for numerical stability.
 *
 *  \param logitp  log(p/(1-p)), the natural (logit) parameter
 *  \param nu   dispersion parameter (>0)
 *  \param n    number of trials (treated as constant, never differentiated)
 *  \return     log Z
 */
template<class Type>
Type calc_logZ(Type logitp, Type nu, int n) {
  using atomic::tiny_ad::isfinite;
  bool ok = (n >= 0 && isfinite(logitp) && isfinite(nu));
  if (!ok) return NAN;
  using atomic::robust_utils::logspace_add;
  Type logZ = Type(-INFINITY);
  Type log_choose_k = Type(0);  // log C(n, 0) = 0
  for (int k = 0; k <= n; ++k) {
    if (k > 0) {
      // log C(n, k) = log C(n, k-1) + log(n-k+1) - log(k)
      log_choose_k += log((double)(n - k + 1)) - log((double)k);
    }
    Type log_term = nu * log_choose_k + (double)k * logitp;
    logZ = logspace_add(logZ, log_term);
  }
  return logZ;
}

/** \brief Conway-Maxwell-Binomial. Calculate mean from natural parameter. */
template<class Type>
Type calc_mean(Type logitp, Type nu, int n) {
  typedef atomic::tiny_ad::variable<1, 1, Type> ADType;
  ADType logitp_(logitp, 0);
  ADType ans = calc_logZ<ADType>(logitp_, nu, n);
  return ans.getDeriv()[0];
}

/** \brief Conway-Maxwell-Binomial. Calculate logit(p) from log(mean).
 *
 *  Inverts E[Y | n, logitp, nu] = exp(mean) for logitp using safeguarded
 *  Newton iteration with bisection backstop. Pure Newton oscillates at
 *  strong overdispersion (nu << 1); the bracket-based variant falls back
 *  to bisection when the Newton step lands outside [logitp_lo, logitp_hi].
 */
template<class Type>
Type calc_logitp(Type mean, Type nu, int n) {
  using atomic::tiny_ad::isfinite;
  bool ok = (n >= 0 && isfinite(mean) && isfinite(nu));
  if (!ok) return NAN;
  int iter_max = 200;
  int expand_max = 200;
  double reltol = 1e-12;
  double abstol = 1e-14;
  typedef atomic::tiny_ad::variable<1, 1, Type> ADType;
  Type mu = mean;
  /* Scaled starting value. At nu = 1 the CMB is binomial and the root is
     exactly logit(mu/n). For nu > 1 it scales roughly as nu * logit(mu/n);
     for nu < 0 it compresses by roughly n. */
  Type logit0 = log(mu / (Type(n) - mu));
  Type x0;
  if (nu >= Type(1)) {
    x0 = nu * logit0;
  } else if (nu > Type(0)) {
    x0 = logit0;
  } else {
    x0 = logit0 / Type(n);
  }
  /* Establish the bracket by outward expansion and verify it before
     iterating. E[Y | logitp] is strictly increasing with limits 0 and n,
     so a root exists for every 0 < mu < n and the expansion terminates. */
  Type logitp_lo = x0;
  Type logitp_hi = x0;
  Type h = fabs(x0);
  if (h < Type(1)) h = Type(1);
  Type m = calc_mean<Type>(x0, nu, n);
  int e = 0;
  bool bracketed = true;
  if (m > mu) {
    while (m > mu && e < expand_max) {
      logitp_lo -= h;
      h *= Type(2);
      e++;
      m = calc_mean<Type>(logitp_lo, nu, n);
    }
    bracketed = (m <= mu);
  } else if (m < mu) {
    while (m < mu && e < expand_max) {
      logitp_hi += h;
      h *= Type(2);
      e++;
      m = calc_mean<Type>(logitp_hi, nu, n);
    }
    bracketed = (m >= mu);
  }
  if (!bracketed) {
    Rf_warning("combinom_utils::calc_logitp: failed to bracket the root");
    return NAN;
  }
  ADType x(x0, 0);
  int i;
  for (i = 0; i < iter_max; i++) {
    x.deriv[0] = Type(1.0);
    ADType y = calc_mean<ADType>(x, nu, n);
    Type residual = y.value - mu;
    /* Apply one Newton correction at convergence. Numerically negligible
       there, but it propagates the implicit derivative of the root: if the
       starting value is already the root (nu == 1), the raw iterate carries
       d(x0)/d(nu) instead of -Cov(Y, lchoose(n,Y))/Var(Y). */
    if (fabs(residual) <= reltol * fabs(mu) || fabs(residual) <= abstol) {
      if (y.deriv[0] > Type(1e-300)) {
        x.value = x.value - residual / y.deriv[0];
      }
      break;
    }
    if (residual > Type(0)) {
      logitp_hi = x.value;
    } else {
      logitp_lo = x.value;
    }
    Type step;
    if (y.deriv[0] > Type(1e-300)) {
      step = -residual / y.deriv[0];
    } else {
      step = (logitp_lo + logitp_hi) / Type(2) - x.value;
    }
    Type x_next = x.value + step;
    if (x_next > logitp_lo && x_next < logitp_hi) {
      x.value = x_next;
    } else {
      x.value = (logitp_lo + logitp_hi) / Type(2);
    }
  }
  if (i == iter_max) {
    Rf_warning("combinom_utils::calc_logitp: Maximum number of iterations exceeded");
  }
  return x.value;
}

} // namespace combinom_utils
