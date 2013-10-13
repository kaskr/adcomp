// Binomial-Poisson mixture model (Royle Biometrics 2004)
// Site-level covariate of p
// Random group effect

DATA_SECTION

  init_int R              // Number of sites
  init_int T              // Number of occasions
  init_int S              // Number of possible values of N
  init_int nG             // Number of groups

  init_ivector N(1,S)         // Possible values of N
  init_ivector nID(1,R)       // Number of observers present at a site
  init_imatrix ID(1,R,1,nID)  // IDs of observer present at each site
  init_imatrix IDind(1,R,1,T) // Group ID RT matrix
  init_vector x(1,R)          // Site-specific covariate
  init_imatrix y(1,R,1,T)     // Count data with R rows and T columns

PARAMETER_SECTION

  init_bounded_number log_lambda(-20.0,20.0,1)
  init_bounded_number p0(-20.0,20.0,1)
  init_bounded_number p1(-20.0,20.0,1)

  init_bounded_number log_sigma(-3.0,3.0,2) // log of random effect SD

  random_effects_vector u(1,nG,2)
  objective_function_value nll

PROCEDURE_SECTION

  for(int i=1;i<=nG;i++)
    prior_N01(u(i));

  for(int i=1;i<=R;i++)
    nll_group(i, p0,p1,log_lambda,log_sigma,u(ID(i)));


SEPARABLE_FUNCTION void prior_N01(const dvariable& u)
   nll += 0.5*square(u);

SEPARABLE_FUNCTION void nll_group(int i, const dvariable& p0,const dvariable& p1,const dvariable& log_lambda, const dvariable& log_sigma, const dvar_vector& u)

  dvariable sigma = exp(log_sigma);
  dvariable lambda = exp(log_lambda);

  double e=1e-12;

  dvar_vector p(1,nID(i));
  dvar_vector f(1,S);
  dvar_vector g(1,S);

  p = 1.0/(1.0+exp(-1.0*(p0 + p1*x(i) + sigma*u)));

  for(int k=1;k<=S;k++)
    f(k) = pow(lambda, N(k)) / exp(lambda + gammln(N(k)+1));

  for(int k=1;k<=S;k++) {
    g(k) = 1.0;
    for(int j=1;j<=T;j++) {
      if(N(k)>=y(i,j))
        g(k) *= exp(gammln(N(k)+1) - gammln(y(i,j)+1) -
                    gammln(N(k)-y(i,j)+1)) * pow(p(IDind(i,j))+e, y(i,j)) *
                    pow(1.0 + e - p(IDind(i,j)), N(k)-y(i,j));
      else
        g(k) = 0.0;
    }
  }

  nll -= log(e + f*g);

TOP_OF_MAIN_SECTION
  arrmblsize = 40000000L;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(30000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(20000000);

