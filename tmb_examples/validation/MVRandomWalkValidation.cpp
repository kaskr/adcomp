// Estimate and validate a multivariate random walk model with correlated increments and correlated observations.
//
// Compare Thygesen et al (submitted, 2016): Validation of state space
// models fitted as mixed effects models
// Casper W. Berg and Kasper Kristensen, 2016

#include <TMB.hpp>

/* Parameter transform */
template <class Type>
Type f(Type x)
{
  return Type(2) / (Type(1) + exp(-Type(2) * x)) - Type(1);
}

template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_ARRAY(obs); /* timeSteps x stateDim */
  PARAMETER(transf_rho);
  PARAMETER(transf_rhoObs);
  PARAMETER_VECTOR(logsds);
  PARAMETER_VECTOR(logsdObs);
  PARAMETER_ARRAY(u); /* State */
  DATA_ARRAY_INDICATOR(keep, obs);
  int timeSteps = obs.dim[1];

  int stateDim = obs.dim[0];
  Type rho = f(transf_rho);
  Type rhoObs = f(transf_rhoObs);
  vector<Type> sds = exp(logsds);
  vector<Type> sdObs = exp(logsdObs);
  // Setup object for evaluating multivariate normal likelihood
  matrix<Type> cov(stateDim, stateDim);
  matrix<Type> covObs(stateDim, stateDim);
  for (int i = 0; i < stateDim; i++)
    for (int j = 0; j < stateDim; j++) {
      cov(i, j) = pow(rho, Type(abs(i - j))) * sds[i] * sds[j];
      covObs(i, j) = pow(rhoObs, Type(abs(i - j))) * sdObs[i] * sdObs[j];
    }
  using namespace density;
  MVNORM_t<Type> neg_log_density(cov);
  MVNORM_t<Type> neg_log_densityObs(covObs);

  /* Define likelihood */
  Type ans = 0;

  // Initial condition
  Type huge = 10;
  for (int i = 0; i < stateDim; i++) ans -= dnorm(u(i, 0), Type(0), huge, true);

  // residuals.setZero();
  for (int i = 1; i < timeSteps; i++)
    ans += neg_log_density(u.col(i) - u.col(i - 1));  // Process likelihood
  for (int i = 0; i < timeSteps; i++) {
    ans += neg_log_densityObs(obs.col(i) - u.col(i),
                              keep.col(i));  // Data likelihood
  }
  ADREPORT(rho);
  ADREPORT(rhoObs);

  return ans;
}
