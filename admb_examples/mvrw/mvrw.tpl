GLOBALS_SECTION 
  #include <df1b2fun.h>
  #include "nLogNormal.h"
DATA_SECTION
  init_int N
  init_int stateDim
  init_matrix obss(1,N,1,stateDim)
PARAMETER_SECTION
  objective_function_value jnll; 
  init_bounded_number rho(0.001,0.999,1);
  init_vector logSdObs(1,stateDim);
  init_vector logSd(1,stateDim);
  random_effects_vector U(1,stateDim*N);
PROCEDURE_SECTION
  
  for(int t=1; t<=(N-1); t++)
     step(t,U((t-1)*stateDim+1,t*stateDim),U(t*stateDim+1,(t+1)*stateDim),logSd,rho);

  for(int t=1; t<=(N-1); t++)
    obs(t,U((t-1)*stateDim+1,t*stateDim),logSdObs);

SEPARABLE_FUNCTION void step(const int t, const dvar_vector& u1,const dvar_vector& u2, const dvar_vector& logSd, const dvariable& rho)
  dvar_matrix fvar(1,stateDim,1,stateDim);
  dvar_matrix fcor(1,stateDim,1,stateDim);
  dvar_vector fsd(1,stateDim);

  fvar.initialize();
  fsd = exp(logSd);
  
  dvar_vector a=u1.shift(1);
  dvar_vector b=u2.shift(1);

  for(int i=1; i<=stateDim; ++i){
        for(int j=1; j<=stateDim; ++j){
          if(i!=j){fcor(i,j)=pow(rho,abs(i-j));}else{fcor(i,j)=1.0;}
        }
  }
  //cout << t << "flaf" << endl;
  fvar=elem_prod(outer_prod(fsd,fsd),fcor);
  jnll+=nLogNormal(a,b,fvar);


SEPARABLE_FUNCTION void obs(const int t, const dvar_vector& u, const dvar_vector& logSdObs)
  dvar_vector var = exp(2.0*logSdObs);
  dvar_vector pred = u.shift(1);		   
  for(int i=1; i<=stateDim; i++){
      //cout << i << endl;
      jnll+=0.5*(log(2.0*M_PI*var(i))+square(obss(t,i)-pred(i))/var(i));
  } 

TOP_OF_MAIN_SECTION
  arrmblsize=2000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(150000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(800000);
  gradient_structure::set_MAX_NVAR_OFFSET(100000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
