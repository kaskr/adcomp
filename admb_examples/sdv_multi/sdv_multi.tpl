// Multitivariat, men strippet ned slik jun vil

DATA_SECTION

  init_int n						// Length of time series
  init_int p						// Dimension of vector

  init_matrix y(1,n,1,p)				// Time series

  int p_n					
  !!p_n = p*n;

  int pp2					
  !!pp2 = p*(p-1)/2;

  int p21					
  !!p21 = p*p-p;

  // Activity files
  !!USER_CODE ad_comm::change_datafile_name("sdv_multi.pha");
  init_int phi_pha
  init_int log_sigma_pha
  init_int mu_x_pha
  init_int off_diag_x_pha
  init_int h_pha


PARAMETER_SECTION

  init_bounded_vector phi(1,p,-.999,.999,phi_pha)				// Correlation coefficient
  init_bounded_vector log_sigma(1,p,-3.0,3.0,log_sigma_pha)		// 0.5*log(variance component)
  init_bounded_vector mu_x(1,p,-10.,3.,mu_x_pha)				// Regression parameters 
  init_bounded_vector off_diag_x(1,pp2,-5.0,3.0,2)		// Off-diagonal in chol(cor(diff))

  sdreport_vector Sigma(1,p)
  sdreport_matrix Sigma_x(1,p,1,p)

  random_effects_vector h(1,p_n,h_pha)					// State variable

  objective_function_value g

PRELIMINARY_CALCS_SECTION
  cout << setprecision(4);

GLOBALS_SECTION

  //#include <df1b2fun.h>
  
PROCEDURE_SECTION

  int i, j;	

  for (j=1;j<=p;j++)
    sf1(log_sigma(j),phi(j),h(j));

  for (i=2;i<=n;i++)
  {
      sf2(log_sigma,phi,h((i-1)*p+1,i*p),h((i-2)*p+1,(i-1)*p));
  }

  for (i=1;i<=n;i++)
  {
    sf3(h((i-1)*p+1,i*p),mu_x,off_diag_x,i);
  }

  if (sd_phase())
  {
    int ii;

    Sigma = exp(2*log_sigma);

    dvar_matrix L2(1,p,1,p);
    L2.initialize();  
    L2(1,1)=1;
    ii=1;
    for (i=1;i<=p;i++)
    {
      L2(i,i)=1;
      for (int j=1;j<i;j++)
      {
        L2(i,j)=off_diag_x(ii++);
      }
      L2(i)(1,i)/=norm(L2(i)(1,i));
    }
    Sigma_x=L2*trans(L2);

  // To avoid numerical instabilities (trying to calculate SD of a constants)
  for(i=1;i<=p;i++)
    Sigma_x(i,i) = exp(log_sigma(i));

  }


SEPARABLE_FUNCTION void sf1(const dvariable& ls,const dvariable& bb,const dvariable& x_1)
  g -= -0.91893853 -ls + 0.5*log(1-square(bb)) - 0.5*square(x_1/mfexp(ls))*(1-square(bb));

SEPARABLE_FUNCTION void sf2(const dvar_vector& ls,const dvar_vector& phi,const dvar_vector& _h_i,const dvar_vector& _h_i1)
  ADUNCONST(dvar_vector,h_i)
  ADUNCONST(dvar_vector,h_i1)
  h_i.shift(1);
  h_i1.shift(1);

  int i, j;

  dvar_vector diff(1,p);
  diff = h_i - elem_prod(phi,h_i1);

  g -= -p*0.91893853 - sum(ls) - .5*norm2(elem_div(diff,mfexp(ls)));

SEPARABLE_FUNCTION void sf3(const dvar_vector& _h_i , const dvar_vector& mu_x, const dvar_vector& Off_diag_x, const int ii)
  ADUNCONST(dvar_vector,h_i)
  h_i.shift(1);

  int i, j;

  dvar_vector log_sigma_y(1,p);
  dvar_vector sigma_y(1,p);
  log_sigma_y = 0.5*(mu_x + h_i);
  sigma_y = mfexp(log_sigma_y);

  dvar_matrix L(1,p,1,p);
  L.initialize();
  L(1,1) = 1.0;
  int k = 1;
  for(i=2;i<=p;i++)
  {
    L(i,i) = 1.0;
    for(j=1;j<=i-1;j++)
      L(i,j) = Off_diag_x(k++);
    L(i)(1,i) /= norm(L(i)(1,i));
  }

  dvariable logdetL = 0.0;
  for(i=1;i<=p;i++)
   logdetL += log(L(i,i));

  g -= -p*0.91893853 - sum(log_sigma_y) - logdetL
	//-.5*norm2(elem_div(solve(L,y(ii)),sigma_y));
	-.5*norm2(solve(L,elem_div(y(ii),sigma_y)));


REPORT_SECTION

TOP_OF_MAIN_SECTION
  arrmblsize = 4000000L;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(500000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2000000);
  gradient_structure::set_MAX_NVAR_OFFSET(50850);


