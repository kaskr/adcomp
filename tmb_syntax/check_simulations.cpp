#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRING(distr);
  DATA_INTEGER(n);
  Type ans = 0;

  if (distr == "norm") {
    PARAMETER(mu);
    PARAMETER(sd);
    vector<Type> x = rnorm(n, mu, sd);
    ans -= dnorm(x, mu, sd, true).sum();
  }
  else if (distr == "gamma") {
    PARAMETER(shape);
    PARAMETER(scale);
    vector<Type> x = rgamma(n, shape, scale);
    ans -= dgamma(x, shape, scale, true).sum();
  }
  else if (distr == "pois") {
    PARAMETER(lambda);
    vector<Type> x = rpois(n, lambda);
    ans -= dpois(x, lambda, true).sum();
  }
  else if (distr == "nbinom") {
    PARAMETER(size);
    PARAMETER(prob);
    vector<Type> x = rnbinom(n, size, prob);
    ans -= dnbinom(x, size, prob, true).sum();
  }
  else if (distr == "nbinom2") {
    PARAMETER(mu);
    PARAMETER(var);
    vector<Type> x = rnbinom2(n, mu, var);
    ans -= dnbinom2(x, mu, var, true).sum();
  }
  else if (distr == "exp") {
    PARAMETER(rate);
    vector<Type> x = rexp(n, rate);
    ans -= dexp(x, rate, true).sum();
  }
  else if (distr == "AR1") {
    PARAMETER(phi);
    vector<Type> x(n);
    density::AR1(phi).simulate(x);
    ans += density::AR1(phi)(x);
  }
  else if (distr == "ARk") {
    PARAMETER_VECTOR(phi);
    vector<Type> x(n);
    density::ARk(phi).simulate(x);
    ans += density::ARk(phi)(x);
  }
  else if (distr == "MVNORM") {
    PARAMETER(phi);
    matrix<Type> Sigma(5,5);
    for(int i=0; i<Sigma.rows(); i++)
      for(int j=0; j<Sigma.rows(); j++)
        Sigma(i,j) = exp( -phi * abs(i - j) );
    density::MVNORM_t<Type> nldens = density::MVNORM(Sigma);
    for(int i = 0; i<n; i++) {
      vector<Type> x = nldens.simulate();
      ans += nldens(x);
    }
  }
  else if (distr == "SEPARABLE") {
    PARAMETER(phi1);
    PARAMETER_VECTOR(phi2);
    array<Type> x(100, 200);
    SEPARABLE( density::ARk(phi2), density::AR1(phi1) ).simulate(x);
    ans += SEPARABLE( density::ARk(phi2), density::AR1(phi1) )(x);
  }
  else error( ("Invalid distribution '" + distr + "'").c_str() );
  return ans;
}
