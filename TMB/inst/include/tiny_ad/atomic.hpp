/********************************************************************
 * Special math functions are available in TMB by including this file
 * in the cpp file.
 ********************************************************************/

/** \file
    \brief Special functions depending on 'tiny_ad'.
*/

template<class Type> Type lgamma(Type x);

namespace atomic {

#define TINY_AD_USE_TINY_VEC 1
#include "tiny_ad/tiny_ad.hpp"
#include "mask.hpp"

/********************************************************************
 * Adding 'pbeta'
 ********************************************************************/
#include "beta/pbeta.hpp"    // Get namespace 'toms708'
TMB_BIND_ATOMIC(pbeta,
		111,
		toms708::pbeta(x[0], x[1], x[2], 1, 0) )

/********************************************************************
 * Adding 'bessel_k'
 ********************************************************************/
#include "bessel/bessel.hpp" // Get namespace 'bessel_utils'
TMB_BIND_ATOMIC(bessel_k,
		11,
		bessel_utils::bessel_k(x[0], x[1], 1.) )

/********************************************************************
 * Adding 'bessel_i'
 ********************************************************************/
#include "bessel/bessel.hpp" // Get namespace 'bessel_utils'
TMB_BIND_ATOMIC(bessel_i,
		11,
		bessel_utils::bessel_i(x[0], x[1], 1.) )

/********************************************************************
 * Adding 'dtweedie'
 ********************************************************************/
#include "tweedie/tweedie.hpp"
TMB_BIND_ATOMIC(log_dtweedie,
		0111,
		tweedie_utils::dtweedie(x[0], x[1], x[2], x[3], true) )

/********************************************************************
 * Adding 'qbeta'
 ********************************************************************/
extern "C" double Rf_qbeta(double, double, double, int, int);
template <class Type>
Type dbeta(Type x, Type shape1, Type shape2) {
  Type logres =
    lgamma(shape1 + shape2) - lgamma(shape1) - lgamma(shape2) +
    (shape1 - 1.) * log(x) + (shape2 - 1.) * log(1. - x);
  return exp(logres);
}
TMB_ATOMIC_VECTOR_FUNCTION(
			   // ATOMIC_NAME
			   qbeta
			   ,
			   // OUTPUT_DIM
			   1
			   ,
			   // ATOMIC_DOUBLE
			   ty[0] = Rf_qbeta(tx[0], tx[1], tx[2], 1, 0);
			   ,
			   // ATOMIC_REVERSE
			   Type p = ty[0];
			   Type a = tx[1];
			   Type b = tx[2];
			   Type tmp = atomic::dbeta(p, a, b);
			   px[0] = 1.0 / tmp * py[0];
			   CppAD::vector<Type> arg(4);
			   arg[0] = p;
			   arg[1] = a;
			   arg[2] = b;
			   arg[3] = Type(1); // 1st order partials wrt. a and b
			   CppAD::vector<Type> D_shape = pbeta(arg);
			   px[1] = -D_shape[1] / tmp * py[0];
			   px[2] = -D_shape[2] / tmp * py[0];
			   )

} // End namespace atomic


/********************************************************************
 * Interfaces
 ********************************************************************/

/** \brief Distribution function of the beta distribution (following R
    argument convention).
    \note Non-centrality parameter (ncp) not implemented.
    \ingroup R_style_distribution
*/
template<class Type>
Type pbeta(Type q, Type shape1, Type shape2){
  CppAD::vector<Type> tx(4);
  tx[0] = q;
  tx[1] = shape1;
  tx[2] = shape2;
  tx[3] = 0; // order
  Type ans = atomic::pbeta(tx)[0];
  return ans;
}

/** \brief Quantile function of the beta distribution (following R
    argument convention).
    \note Non-centrality parameter (ncp) not implemented.
    \ingroup R_style_distribution
*/
template<class Type>
Type qbeta(Type p, Type shape1, Type shape2){
  CppAD::vector<Type> tx(3);
  tx[0] = p;
  tx[1] = shape1;
  tx[2] = shape2;
  Type ans = atomic::qbeta(tx)[0];
  return ans;
}

/** \brief bessel_k function (same as besselK from R).
    \note Derivatives wrt. both arguments are implemented
    \ingroup special_functions
*/
template<class Type>
Type bessel_k(Type x, Type nu){
  CppAD::vector<Type> tx(3);
  tx[0] = x;
  tx[1] = nu;
  tx[2] = 0;
  Type ans = atomic::bessel_k(tx)[0];
  return ans;
}

/** \brief bessel_i function (same as besselI from R).
    \note Derivatives wrt. both arguments are implemented
    \ingroup special_functions
*/
template<class Type>
Type bessel_i(Type x, Type nu){
  CppAD::vector<Type> tx(3);
  tx[0] = x;
  tx[1] = nu;
  tx[2] = 0;
  Type ans = atomic::bessel_i(tx)[0];
  return ans;
}

/** \brief dtweedie function (same as dtweedie.series from R package
    'tweedie').

    Silently returns NaN if not within the valid parameter range:
    \f[ (0 \leq y) \land (0 < \mu) \land (0 < \phi) \land (1 < p) \land (p < 2) \f] .

    \note Parameter order differs from the R version.

    \warning The derivative wrt. the y argument is disabled
    (zero). Hence the tweedie distribution can only be used for *data*
    (not random effects).

    \ingroup R_style_distribution
*/
template<class Type>
Type dtweedie(Type y, Type mu, Type phi, Type p, int give_log = 0) {
  CppAD::vector<Type> tx(5);
  tx[0] = y;
  tx[1] = mu;
  tx[2] = phi;
  tx[3] = p;
  tx[4] = 0;
  Type ans = atomic::log_dtweedie(tx)[0];
  return ( give_log ? ans : exp(ans) );
}
