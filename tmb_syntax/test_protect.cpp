// Trigger PROTECT bug in TMB version 1.7.8
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  PARAMETER(p);

  if (p>0) {
    int n = 10000;
    matrix<Type> qwerzxcv(n,n);
    qwerzxcv.setIdentity();
    // Explanation:
    // * 'qwerzxcv' symbol not in the symbol table
    // * Matrix is large so install('qwerzxcv') triggers garbage collector
    REPORT(qwerzxcv);
  }

  return 0;
}
