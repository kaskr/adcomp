#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
    DATA_VECTOR(y);
    PARAMETER(mu);
    PARAMETER(logSigma);
    PARAMETER_VECTOR(s); // saddlepoints
    
    // Return K_y(s)-s^T y
    // K_y(s) = sum K_{y[i]}(s[i])
    // K_N(mu, logSigma)(s) = mu*s + sigma^2 * s^2 / 2
    
    Type sigma = exp(logSigma);
    int n = y.size();
    
    // Build CGF
    Type K = 0;
    for(int i=0; i<n; i++){
        K += mu * s(i) + 0.5*s(i)*s(i)*sigma*sigma;
    }
    
    // Build inner problem
    Type res = K - (s*y).sum();
    
    // report sigma
    ADREPORT(sigma);
    
    return res;
}