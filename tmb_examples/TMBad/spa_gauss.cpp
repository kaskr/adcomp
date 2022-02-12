// Demonstrate saddlepoint approximation (SPA)
#include <TMB.hpp>

// Define the ad type needed
typedef TMBad::ad_aug ad;

struct spa_gauss {
  vector<ad> y; // Data
  ad mu, sigma; // Parameters
  // Return K_y(s)-s^T y
  // K_y(s) = sum K_{y[i]}(s[i])
  // K_N(mu, sigma)(s) = mu*s + sigma^2 * s^2 / 2
  ad operator()(vector<ad> s) {
    int n = s.size();
    // Build CGF
    ad K = 0;
    for(int i=0; i<n; i++){
      K += mu * s(i) + 0.5 * s(i) * s(i) * sigma * sigma;
    }
    // Build inner problem
    ad res = K - (s*y).sum();
    return res;
  }
  // Type can by 'ad' or 'double'
  template<class Type>
  Type eval_nldens(vector<Type> &start) {
    newton::newton_config cfg;
    cfg.SPA = true;
    cfg.sparse = true;
    cfg.trace  = false; // Trace inner optimize?
    Type res = newton::Laplace(*this, start, cfg);
    return res;
  }
};

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_VECTOR(s); // saddlepoint initial guess
  PARAMETER(mu);
  PARAMETER(logSigma);
  Type sigma = exp(logSigma);

  spa_gauss spa = {y, mu, sigma};

  // s = start guess on input
  // overwritten with solution on output
  Type res = spa.eval_nldens(s);

  // report saddlepoint and sigma
  REPORT(s); // sd not needed
  ADREPORT(sigma);

  return res;
}
