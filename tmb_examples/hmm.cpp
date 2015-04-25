// Inference in a 'double-well' stochastic differential equation using HMM filter.
#include <TMB.hpp>
#include "hmm_filter.hpp"

/*
  User-specified SDE of the form:

  dX_t = f(X_t) * dt + g(X_t) * dB_t

*/
template<class Type>
struct sde_t{
  Type lambda, gamma, sigmaX;
  sde_t(Type lambda_, Type gamma_, Type sigmaX_){
    lambda = lambda_; gamma = gamma_; sigmaX = sigmaX_;
  }
  // f(x)
  Type advection  (Type x){ return lambda * x - gamma * pow(x,3); }
  // g(x)
  Type dispersion (Type x){ return sigmaX; }
};


template<class Type>
Type objective_function<Type>::operator() ()
{
  // Grid of the state space
  DATA_VECTOR(grid);
  // Measurements are taken at equidistant times
  DATA_SCALAR(dt);
  // Measurements are grid-cell pointers (zero-based)
  DATA_IVECTOR(yobs);

  PARAMETER(lambda);
  PARAMETER(gamma);
  PARAMETER(logsX);
  PARAMETER(logsY);
  Type sigmaX=exp(logsX);
  Type sigmaY=exp(logsY);

  /* Construct the SDE */
  sde_t<Type> sde(lambda, gamma, sigmaX);

  /* Finite volume disretize and report generator */
  fvade_t<sde_t, Type> fvol = fvade(sde, grid);
  REPORT(fvol.A);

  /* Construct likelihood function */
  hmm_filter<Type> hmm_nll(fvol, dt);
  hmm_nll.setGaussianError(sigmaY);
  Type ans = hmm_nll(yobs);

  return ans;
}

