// Demonstrate adaptive solver of TMBad
#include <TMB.hpp>

/*
  This class holds the function we wish to minimize

  - It holds a mix of data and parameters
  - We don't want to specify which is which

*/
template<class Type>
struct Functor {
  matrix<Type> m;
  vector<Type> v;
  Functor(const matrix<Type> &m, const vector<Type> &v) : m(m), v(v) {}
  Type operator()(const vector<Type> &x) {
    v = v * v;
    vector<Type> y = x - v;
    vector<Type> my = m * y;
    return (y * my).sum() + exp(x).sum();
  }
};

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace newton;
  DATA_MATRIX(m);
  PARAMETER_VECTOR(x);
  DATA_STRUCT(cfg, newton_config_t);
  Functor<TMBad::ad_aug> F(m, x);
  vector<Type> sol = Newton(F, x, cfg);
  ADREPORT(sol);
  return sol.sum();
}
