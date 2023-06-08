#include <unsupported/Eigen/FFT>

namespace atomic {

template<bool adjoint=false>
void fft_work(const CppAD::vector<double> &x,
              CppAD::vector<double> &y) {
  int ncplx = x.size() / 2;
  typedef std::complex<double> C;
  C* X = (C*) x.data();
  C* Y = (C*) y.data();
  Eigen::FFT<double> f;
  // How Eigen/FFT differs: invertible scaling is performed
  f.SetFlag(f.Unscaled);
  if (!adjoint)
    f.fwd(Y, X, ncplx);
  else
    f.inv(Y, X, ncplx);
}

TMB_ATOMIC_VECTOR_FUNCTION_DECLARE(fft)  // So can be used by ifft
TMB_ATOMIC_VECTOR_FUNCTION_DECLARE(ifft) // So can be used by fft
TMB_ATOMIC_VECTOR_FUNCTION_DEFINE( fft, tx.size(), fft_work<0>(tx, ty), px = ifft(py))
TMB_ATOMIC_VECTOR_FUNCTION_DEFINE(ifft, tx.size(), fft_work<1>(tx, ty), px =  fft(py))

/** \brief FFT (unscaled)
    \param xc Complex vector to be transformed
    \param inverse Apply the (unscaled) inverse transform?
    This interface is identical to 'stats::fft' so does **not** perform inverse scaling.
    \ingroup matrix_functions
    \note To use this function the header must be manually included.
    \note The default implementation is not the most efficient available. The more efficient FFTW library can be used by setting the preprocessor flag `EIGEN_FFTW_DEFAULT`.
*/
template<class Type>
vector<std::complex<Type> > fft(vector<std::complex<Type> > xc,
                                bool inverse = false) {
  CppAD::vector<Type> x(2 * xc.size());
  Type* px = x.data();
  Type* pxc = (Type*) xc.data();
  for (size_t i=0; i<x.size(); i++) {
    px[i] = pxc[i];
  }
  CppAD::vector<Type> y;
  if (!inverse)
    y = atomic::fft(x);
  else
    y = atomic::ifft(x);
  px = (Type*) y.data();
  for (size_t i=0; i<x.size(); i++) {
    pxc[i] = px[i];
  }
  return xc;
}

}
