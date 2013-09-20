// Ordinal response model with random effects. 
// The SOCATT data used in the software comparison:
// http://www.cmm.bristol.ac.uk/learning-training/multilevel-m-software/index.shtml


DATA_SECTION

  init_int n			// Number of observations
  init_ivector y(1,n)		// Response vector
  init_int S			// Number of response categories
  init_int p			// Number of fixed effects
  init_matrix X(1,n,1,p)	// Fixed effects design matrix
  init_int M			// Number of random effects
  init_ivector ngroup(1,M)      // Group indicator

PARAMETER_SECTION

  init_bounded_vector b(1,p,-5,5) 		// Fixed effects
  init_bounded_number logsigma(-3.0,2.0)	// 0.5*log-variance component
  init_bounded_vector tmpk(1,S-1,-6,6)		// Variables underlying kappa

  random_effects_vector u(1,M)			// Random effects

  objective_function_value g

PROCEDURE_SECTION

  int ii = 1;

  // Main loop calculating likelihood contribution from each cluster
  for(int i=1;i<=(int) M;i++)
    sep_fun(i,ii,u(i),b,logsigma,tmpk);


SEPARABLE_FUNCTION void sep_fun(int& i, int& ii, const dvariable& u,const dvar_vector& b,const dvariable& logsigma,const dvar_vector& tmpk)

  dvariable sigma = exp(logsigma);

  dvar_vector alpha(1,S-1);
  alpha(1) = tmpk(1);
  for(int s=2;s<=S-1;s++)
    alpha(s) = alpha(s-1) + mfexp(tmpk(s));

  dvariable P, eta;                      

  g -= -0.5*log(2*3.1415927) - 0.5*square(u);           // 

  for(int j=1;j<=ngroup(i);j++)
  {
    eta = X(ii)*b + sigma*u;             		// X(i) is i'th row of X
    if(y(ii)==S)				
      P = 1.0;
    else
      P = 1/(1+exp(-(alpha(y(ii))-eta)));
    if(y(ii) > 1)
      P -= 1/(1+exp(-(alpha(y(ii)-1)-eta)));
    g -= log(1.e-20+P);                   		

    ii++;
  }

TOP_OF_MAIN_SECTION
  arrmblsize = 4000000L;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(200000);
  gradient_structure::set_MAX_NVAR_OFFSET(2000);

