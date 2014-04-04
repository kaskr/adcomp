DATA_SECTION
  init_int nobs
  init_matrix Yx(1,nobs,1,2)
  vector Y(1,nobs);
  !! Y=column(Yx,1);
  vector x(1,nobs);
  !! x=column(Yx,2);

PARAMETER_SECTION
  init_number a
  init_number b
  init_number logSigma
  sdreport_number sigmasq
  objective_function_value nll

PROCEDURE_SECTION
  sigmasq=exp(2*logSigma);
  nll=0.5*(nobs*log(2*M_PI*sigmasq)+sum(square(Y-(a+b*x)))/sigmasq);

TOP_OF_MAIN_SECTION
  arrmblsize=50000000;
 