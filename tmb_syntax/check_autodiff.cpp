// Check correctness of 'autdiff' namespace
#include <TMB.hpp>

// Function with one output
struct func1 {
  template <class T>
  T operator() (vector<T> x) {  // Evaluate function
    return x.prod();
  }
};

// Function with three output
struct func3 {
  template <class T>
  vector<T> operator()(vector<T> x) {  // Evaluate function
    vector<T> y(3);
    y(0) = x(0);    
    y(1) = x.sum();
    y(2) = x.prod();
    return y;
  }
};

// Function with local variables
// grad_x(func2) = [ exp(a), exp(b) ]
// jac_ab = diag( [ exp(a), exp(b) ] )
template<class Type>
struct func2 {
  Type a;
  Type b;
  template <class T>
  T operator()(vector<T> x) {  // Evaluate function
    T c = exp(a) * x[0] + exp(b) * x[1];
    return c;
  }
};

template<class Type>
Type objective_function<Type>::operator() () {
  PARAMETER_VECTOR(theta);
  DATA_INTEGER(select);

  if (select == 1) {
    func1 f;
    // Calculate gradient and hessian
    vector<Type> g = autodiff::gradient(f, theta);
    matrix<Type> h = autodiff::hessian (f, theta);
    REPORT(g);
    REPORT(h);
    ADREPORT(f(theta));
  }

  if (select == 2) {
    func2<Type> f = { theta[0], theta[1] };
    // Calculate gradient and hessian
    vector<Type> g = autodiff::gradient(f, theta);
    ADREPORT(g);
  }

  if (select == 3) {  
    func3 f;
    // Calculate jacobian
    matrix<Type> j = autodiff::jacobian(f, theta);
    REPORT(j);
    ADREPORT(f(theta));
  }

  // Exit
  return 0;
}
