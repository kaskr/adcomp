// Adaptive integration using 'tiny_ad'
#include <TMB.hpp>

namespace my_atomic {
  /*
    This class evaluates the marginal density of x where
    u     ~ Normal( mu, sd^2 )
    x | u ~ Binom ( n , plogis(u) )
  */
  template<class Float>
  struct GaussBinomial_t {
    typedef Float Scalar; // Required by integrate
    Float x, n;           // Data
    Float mu, sd;         // Parameters
    // Evaluate joint density of (u, x)
    Float operator() (Float u) {
      Float ans = 0;
      ans += dnorm(u, Float(0.), Float(1.), true);
      Float p = invlogit(sd * u + mu);
      p = squeeze(p);
      ans += x * log(p) + (n - x) * log(1. - p);
      ans = exp(ans);
      // Avoid NaNs in the gradient:
      if (ans == 0) ans = 0;
      // Avoid NaNs in the tail of the distribution:
      using atomic::tiny_ad::isfinite;
      if (!isfinite(ans)) ans = 0;
      return ans;
    }
    // Integrate latent variable (u) out
    Float marginal() {
      using gauss_kronrod::integrate;
      Float ans =
	integrate(*this, -INFINITY, INFINITY);
      if(n > 1) {
	using atomic::gamma_utils::gammafn;
	ans /= gammafn(n+1) / (gammafn(x+1) * gammafn(n-x+1));
      }
      return ans;
    }
  };

  // ****** How to use it in TMB:
  // 1. Create an evaluator 'eval' for previous class
  template<class Float>
  Float eval(Float x, Float n, Float mu, Float sd) {
    GaussBinomial_t<Float> f = {x, n, mu, sd};
    return f.marginal();
  }
  // 2. Run 'eval' through tiny_ad and obtain an atomic function
  //    'func'.  The '0011' tells tiny_ad that we only need
  //    derivatives wrt. mu and sd.
  TMB_BIND_ATOMIC(func, 0011, eval(x[0], x[1], x[2], x[3]))
  // 3. Create a more user-friendly version ('func' takes vector
  //    arguments and there's a final invisible argument that
  //    corresponds to the derivative order)
  template<class Type>
  Type GaussBinomial(Type x, Type n, Type mu, Type sd) {
    vector<Type> args(5); // Last index reserved for derivative order
    args << x, n, mu, sd, 0;
    return my_atomic::func(CppAD::vector<Type>(args))[0];
  }
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);
  DATA_VECTOR(n);
  DATA_MATRIX(A);
  PARAMETER_VECTOR(b);
  vector<Type> mu = A * b;
  PARAMETER(logsd);
  Type sd = exp(logsd);
  Type ans = 0;
  Type tiny = 0.0; // Set to 1e-12 for robustness
  for(int i=0; i < x.size(); i++) {
    ans -= log( my_atomic::GaussBinomial(x(i), n(i), mu(i), sd) + tiny );
  }
  return ans;
}
