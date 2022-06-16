// Demonstrate 2D interpolation operator
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using tmbutils::interpol2D;
  using tmbutils::interpol2D_config;
  // Data for interpolation
  DATA_MATRIX(A);
  DATA_VECTOR(x_range);
  DATA_VECTOR(y_range);
  // Test evaluation for theese values
  PARAMETER_VECTOR(x);
  PARAMETER_VECTOR(y);  
  DATA_SCALAR(R);
  interpol2D_config<Type> cfg(R);
  interpol2D<Type> op(A, x_range, y_range, cfg);
  vector<Type> f = op(x, y);
  REPORT(f);
  return f.sum();
}
