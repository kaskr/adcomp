#include <TMB.hpp>
#include <fft.hpp> // atomic fft

template<class T>
vector<std::complex<T> > real2cplx(vector<T> x) {
  vector<std::complex<T> > xc(x.size());
  for (int i=0; i<x.size(); i++) xc[i] = x[i];
  return xc;
}
template<class T>
vector<T> cplx2real(vector<std::complex<T> > xc) {
  vector<T> x(xc.size());
  for (int i=0; i<x.size(); i++) x[i] = xc[i].real();
  return x;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // y ~ MVNORM(0, Sigma)
  PARAMETER(rho);
  DATA_VECTOR(x);
  DATA_VECTOR(y);
  vector<Type> C = exp(-rho*x);
  // fft(y) ~ N(0, F*Sigma*F^H)
  vector<std::complex<Type> > C2 = atomic::fft(real2cplx(C));
  vector<std::complex<Type> > sd = sqrt(C2);
  vector<std::complex<Type> > y2 = atomic::fft(real2cplx(y));
  Type f = -dnorm(y2, std::complex<Type>(0), sd, true).sum().real();
  return f;
}
