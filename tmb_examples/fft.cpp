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
  vector<Type> sd = atomic::fft(cplx(C)).real().sqrt();
  vector<std::complex<Type> > Fx = atomic::fft(cplx(x));
  vector<Type> y = (Fx * Fx.conjugate()).real().sqrt(); // modulus
  y = y / sqrt((Type) y.size());
  return dnorm(y, Type(0), sd, true).sum();
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
