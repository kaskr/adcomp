DATA_SECTION
  init_int N
  init_vector Y(1,N)

PARAMETER_SECTION
  init_number logr0
  init_number logtheta
  init_bounded_number logK(4.6,7.6)
  init_number logQ
  init_number logR

  random_effects_vector X(1,N);

  sdreport_number r0
  sdreport_number theta
  sdreport_number K
  sdreport_number Q
  sdreport_number R

  objective_function_value jnll
  			   
PROCEDURE_SECTION
  for(int i=2; i<=N; ++i){
    step(X(i-1),X(i),logr0,logK,logtheta,logQ);  }

  for(int i=1; i<=N; ++i){ 
    obs(X(i),logR,i);}

  r0=exp(logr0);  theta=exp(logtheta); K=exp(logK); Q=exp(logQ); R=exp(logR);

SEPARABLE_FUNCTION void step(const dvariable& x1, const dvariable& x2, const dvariable& logr0, const dvariable& logK, const dvariable& logtheta, const dvariable& logQ)
  dvariable var=exp(logQ);
  dvariable m=x1 + exp(logr0) * (1.0 - pow(exp(x1)/exp(logK),exp(logtheta))); 
  jnll+=0.5*(log(2.0*M_PI*var)+square(x2-m)/var);

SEPARABLE_FUNCTION void obs(const dvariable& x, const dvariable& logR, int i)
  dvariable var=exp(logR);
  jnll+=0.5*(log(2.0*M_PI*var)+square(x-Y(i))/var);
