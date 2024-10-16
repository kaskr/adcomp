namespace robust_utils {

  using namespace atomic::tiny_ad;

  // logspace_add
  // template<class T>
  // T fmax2(T x, T y) {return (x < y) ? y : x;}
  template<class T>
  T logspace_add (const T &logx, const T &logy) {
    // Was:
    //  fmax2 (logx, logy) + log1p (exp (-fabs (logx - logy)));
    return ( logx < logy ?
             logy + log1p (exp (logx - logy)) :
             logx + log1p (exp (logy - logx)) );
  }
  // logspace_sub
  template<class T>
  T R_Log1_Exp (const T &x) {
    return ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x))) ;
  }
  template<class T>
  T logspace_sub (const T &logx, const T &logy) {
    return logx + R_Log1_Exp(logy - logx);
  }

  /** \brief Robust binomial distribution *without* normalizing constant

      Parameterized via logit(p).

      \internal

      x = logit(p) = log( p / (1-p) )

      p   = 1/(1+exp(-x))
      1-p = 1/(1+exp( x))

      log(p)   = -log( exp(0) + exp(-x) ) = -logspace_add( 0 , -x )
      log(1-p) = -log( exp(0) + exp( x) ) = -logspace_add( 0 ,  x )

  */
  template<class Float>
  Float dbinom_robust(Float k, Float size, Float logit_p, int give_log=0)
  {
    Float zero = 0;
    Float log_p   = -logspace_add( zero , Float(-logit_p) );
    Float log_1mp = -logspace_add( zero , logit_p );
    Float logres = k * log_p + (size-k) * log_1mp;
    return ( give_log ? logres : exp(logres) );
  }

}
