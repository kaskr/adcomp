// Demonstrate how to change data without re-taping
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(a); DATA_UPDATE(a);
  DATA_ARRAY (b); DATA_UPDATE(b);
  DATA_MATRIX(c); DATA_UPDATE(c);
  DATA_SCALAR(d); DATA_UPDATE(d);
  PARAMETER  (x);

  Type sd = a.sum() + b.sum() + c.sum() + d;
  Type nll = -dnorm(x, Type(0), sd, true);
  return nll;
}
