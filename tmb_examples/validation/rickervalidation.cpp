// Estimate and validate a Ricker model based on data simulated from the logistic map
//
// Compare Thygesen et al (submitted, 2016): Validation of state space models
// fitted as mixed effects models
//
// This file implements the "Theta logistic population model" from
// Pedersen et al 2012, Ecol. Modelling. With theta=1, this is the Ricker
// model.
//
// Uffe Høgsbro Thygesen and Kasper Kristensen, 2016

#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  /* Data section */
  DATA_VECTOR(Y);                  // Counted abundance
  DATA_VECTOR_INDICATOR(keep, Y);  // For one-step predictions

  /* Parameter section */
  PARAMETER_VECTOR(X);  // Latent states. As last as long as Y;
                        // extra elements are not used
  PARAMETER(logr);      // Growth rate
  PARAMETER(logtheta);  // With theta=1, the Ricker model
  PARAMETER(logK);      // Carrying capacity
  PARAMETER(logQ);      // Process noise
  PARAMETER(logS);      // Sample size controlling measurement noise

  /* Procedure section */

  Type r = exp(logr);
  Type theta = exp(logtheta);
  Type K = exp(logK);
  Type Q = exp(logQ);
  Type S = exp(logS);

  int timeSteps = Y.size();
  Type nll = 0;

  // Contributions from state transitions
  for (int i = 1; i < timeSteps; i++) {
    Type m = X[i - 1] + r * (1.0 - pow(exp(X[i - 1]) / K, theta));
    nll -= dnorm(X[i], m, sqrt(Q), true);
  }

  // Contributions from observations
  for (int i = 0; i < timeSteps; i++) {
    nll -= keep(i) * dpois(Y[i], S * exp(X[i]),
                           true);  // keep(i) for one-step predictions
  }

  return nll;
}
