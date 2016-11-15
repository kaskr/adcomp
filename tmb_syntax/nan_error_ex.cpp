// Illustrates how to make the debugger catch a floating point error.
#include <TMB.hpp>
#include <fenv.h> // Extra line needed

template<class Type>
Type objective_function<Type>::operator() ()
{
  feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO | FE_UNDERFLOW); // Extra line needed
	
  DATA_SCALAR(lambda);
  PARAMETER(x);
  Type f;
  f = sqrt(-1.);        // FE_INVALID   ( sqrt(-1.) returns NaN )
  //f = 1./0.;          // FE_DIVBYZERO ( division by zero )
  //f = exp(100000.);   // FE_OVERFLOW  ( exp(100000.) returns Inf )   [Does not work on all platforms]
  //f = exp(-100000.);  // FE_UNDERFLOW ( exp(-100000.) returns 0 )
  return f;
}
