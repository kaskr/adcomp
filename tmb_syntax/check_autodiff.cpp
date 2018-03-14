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
