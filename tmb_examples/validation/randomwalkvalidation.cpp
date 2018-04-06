// Estimate and validate a random walk model with and without drift
//
// Compare Thygesen et al (submitted, 2016): Validation of state space models
// fitted as mixed effects models
//
// Uffe Høgsbro Thygesen and Kasper Kristensen, 2016

#include <TMB.hpp>
template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_VECTOR(y);                  // Observations
  DATA_VECTOR_INDICATOR(keep, y);  // For one-step predictions

  DATA_SCALAR(huge);
  PARAMETER_VECTOR(x);
  PARAMETER(mu);
  PARAMETER(logsigma);
  PARAMETER(logs);

  // Initial condition
  Type nll = -dnorm(x(0), Type(0), huge, true);

  // Increments
  for (int i = 1; i < x.size(); ++i)
    nll -= dnorm(x(i), x(i - 1) + mu, exp(logsigma), true);

  // Observations
  for (int i = 0; i < y.size(); ++i)
    nll -= keep(i) * dnorm(y(i), x(i), exp(logs), true);

  return nll;
}
