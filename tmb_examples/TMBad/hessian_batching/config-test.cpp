#include <TMB.hpp>
template <class Type>
Type objective_function<Type>::operator()() {
  PARAMETER(mu);
  PARAMETER_VECTOR(u);
  Type ans = 0, s = 0;
  for (int i = 0; i < u.size(); ++i) {
    s += u[i];
    ans += Type(.5) * u[i] * u[i] + Type(.01) * exp(s + mu);
  }
  return ans + mu * mu;
}
