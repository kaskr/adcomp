#include <TMB.hpp>
template <class Type>
Type objective_function<Type>::operator() ()
{
    DATA_VECTOR(x);
    PARAMETER(k);
    PARAMETER(lambda);
    
    Type jnll, dnchi, u;
    jnll = dnchi = 0;
    u = k / 2 - 1;
    for (int i = 0; i < x.size(); i++) {
        dnchi = 0.5*exp(-(x[i] + lambda) / 2) * pow(Type(x[i] / lambda), Type(u / 2)) *besselI2(sqrt(lambda*x[i]), u);
        jnll -= log(dnchi);
        dnchi = 0;
    }
    return jnll;
}