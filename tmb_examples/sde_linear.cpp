// Inference in a linear scalar stochastic differential equation.
//
// dX = - lambda*X*dt + sigmaX*dB
//
// based on discrete observations
//
// Y(i) = X(t(i)) + e(i)
//
// where e(i) is N(0,sigmaY^2)
//
// Latent variables are the states.
//
// We use Euler approximation to evalaute transition densities. The time mesh for this
// discretization is finer than the sample interval, i.e. some (many) states are unobserved.

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(tsim);     // Time points where X is simulated 
  DATA_VECTOR(iobs);     // Indeces into tsim where X is observed
  DATA_VECTOR(Y);        // Observations taken. Must have same length as iobs.

  PARAMETER_VECTOR(X);   // States at tsim. Length = length(tsim)
  PARAMETER(lambda);     // Rate parameter in the SDE				     
  PARAMETER(logsX);      // log(sigmaX) where sigmaX is noise intensity in the SDE   
  PARAMETER(logsY);      // log(sigmaY) where sigmaY is std.dev. on measurement error

  Type sX=exp(logsX);
  Type sY=exp(logsY);

  Type ans=0;  // ans will be the resulting likelihood

  vector<Type> dt(tsim.size()-1);

  for(int i=0;i<dt.size();i++)   
    dt(i) = tsim(i+1)-tsim(i);

  // Include likelihood contributions from state transitions
  for(int i=0;i<dt.size();i++){
    ans -= dnorm(X(i+1),X(i) + lambda*X(i)*dt(i),sX*sqrt(dt(i)),1);
  }

  // Include likelihood contributions from measurements
  for(int i=0;i<Y.size(); ++i){
    int j=CppAD::Integer(iobs(i));
    ans-=dnorm(Y(i),X(j),sY,1);
  }

  return ans;
}
