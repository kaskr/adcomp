// Similar to example 'adaptive_integration' using CppAD Romberg integration. REGISTER_ATOMIC is used to reduce tape size.
#include <TMB.hpp>

template<class Type>
struct univariate {
  Type x, n, mu, sd; // Parameter in integrand
  Type operator() (Type u) {
    Type ans = 0;
    ans += dnorm(u, Type(0.), Type(1.), true);
    Type p = invlogit(sd * u + mu);
    p = squeeze(p);
    ans += x * log(p) + (n - x) * log(1. - p);
    ans = exp(ans);
    return ans;
  }
};

/*
  This function evaluates the marginal density of x where
  u     ~ Normal( mu, sd^2 )
  x | u ~ Binom ( n , plogis(u) )
*/
template<class Type>
vector<Type> GaussBinomial(vector<Type> input) {
  Type x  = input[0], n  = input[1];         // Data
  Type mu = input[2], sd = input[3];         // Parameters
  univariate<Type> f = {x, n, mu, sd};
  Type a = -5, b = 5;
  vector<Type> res(1);
  res[0] = romberg::integrate(f, a, b);
  res[0] /= exp(lgamma(n+1) - (lgamma(x+1) + lgamma(n-x+1)));
  return res;
}

REGISTER_ATOMIC(GaussBinomial)

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
    vector<Type> input(4);
    input << x(i), n(i), mu(i), sd;
    ans -= log( GaussBinomial(input)[0] + tiny );
  }
  return ans;
}
