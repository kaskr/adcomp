// Multivariate normal distribution with circulant covariance
#include <TMB.hpp>
#include <fft.hpp> // atomic fft

// Helper to convert real vector to complex vector
template<class T>
vector<std::complex<T> > cplx(vector<T> x) {
  vector<std::complex<T> > xc(x.size());
  // Careful to set imag=0. T() is not zero by default!
  for (int i=0; i<x.size(); i++) xc[i] = x[i];
  return xc;
}

// dmvnorm for circulant covariance
template<class Type>
Type log_dmvnorm_fft(vector<Type> x, vector<Type> C) {
  vector<std::complex<Type> > sd = atomic::fft(cplx(C)).sqrt();
  vector<std::complex<Type> > y = atomic::fft(cplx(x), true);
  std::complex<Type> zero(0);
  return
    dnorm(y, zero, sd, true).sum().real() + .5*log((double)C.size());
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // y ~ MVNORM(0, Sigma)
  PARAMETER(rho);
  DATA_VECTOR(d); // circ distance
  DATA_VECTOR(x); // observation
  vector<Type> C = exp(-rho*d);
  return -log_dmvnorm_fft(x, C);
}
