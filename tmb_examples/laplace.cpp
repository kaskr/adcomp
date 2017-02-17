// Laplace approximation from scratch demonstrated on 'spatial' example.
// The purpose is to show how the namespace 'autodiff' can be used to
// differentiate an iterative solver (the inner problem of the Laplace
// approximation). It is not meant to replace the approach used in the
// 'spatial' example (!)

#include <TMB.hpp>

/* Generic Laplace approximation (everywhere convex problems only). */
template<class Type, class Functor>
struct laplace_t {
  Functor f;       // User's implementation of joint nll
  vector<Type>& u; // Random effect vector
  int niter;       // Number of Newton iterations
  laplace_t(Functor f_, vector<Type> &u_, int niter_) :
    f(f_), u(u_), niter(niter_) {}
  Type operator()(){
    // Solve inner problem - Newton iterations
    for (int i=0; i<niter; i++){
      vector<Type> g = autodiff::gradient(f, u);
      matrix<Type> H = autodiff::hessian(f, u);
      u = u - atomic::matinv(H) * g;
    }
    // Laplace approximation
    matrix<Type> H = autodiff::hessian(f, u);
    Type ans = .5 * atomic::logdet(H) + f(u);
    ans -= .5 * Type(u.size()) * log(2.0 * M_PI);
    return ans;
  }
};
template<class Type, class Functor>
Type laplace(Functor f, vector<Type> &u, int niter){
  laplace_t<Type, Functor> L(f, u, niter);
  return L();
}

/* The following is (almost) copy-pasted from the 'spatial' example */
template<class Type>
struct joint_nll {

  /* Data and parameter objects for spatial example: */
  vector<Type> y;
  matrix<Type> X;
  matrix<Type> dd;
  vector<Type> b;
  Type a;
  Type log_sigma;

  /* Constructor */
  joint_nll(vector<Type> y_,
	    matrix<Type> X_,
	    matrix<Type> dd_,
	    vector<Type> b_,
	    Type a_,
	    Type log_sigma_) :
    y(y_), X(X_), dd(dd_), b(b_),
    a(a_), log_sigma(log_sigma_) {}

  /* Evaluate the negative joint log-likelihood as function of the
     random effects */
  template <typename T>
  T operator()(vector<T> u) {
    int n = u.size();
    T res=0;
    vector<T> eta = T(exp(log_sigma)) * u;
    vector<Type> mu = X * b;
    eta = eta + mu.template cast<T>();
    matrix<T> cov(n,n); 
    for (int i=0; i<n; i++)
      {
	cov(i,i) = 1.0;
	for (int j=0; j<i; j++)
	  {
	    // Exponentially decaying correlation
	    cov(i,j) = exp(-a * dd(i,j));
	    cov(j,i) = cov(i,j);
	  }
      }
    density::MVNORM_t<T> neg_log_density(cov);
    res += neg_log_density(u);
    // logdpois = N log lam - lam
    for(int i=0; i<y.size(); i++)
      res -= T(y[i]) * eta[i] - exp(eta[i]);
    return res;
  }
};


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_MATRIX(X)
  DATA_MATRIX(dd)
  PARAMETER_VECTOR(b);
  PARAMETER(a);
  PARAMETER(log_sigma);
  int n = dd.rows();

  // Construct joint negative log-likelihood
  joint_nll<Type> jnll(y, X, dd, b, a, log_sigma);

  // Random effect initial guess
  vector<Type> u(n);
  u.setZero();

  // Calculate Laplace approx (updates u)
  DATA_INTEGER(niter);
  Type res = laplace(jnll, u, niter);
  ADREPORT(u)

  return res;
}
