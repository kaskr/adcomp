// State space assessment model from Nielsen and Berg 2014, Fisheries Research.
//  --------------------------------------------------------------------------
// Copyright (c) 2014, Anders Nielsen <an@aqua.dtu.dk>, 
// Casper Berg <cbe@aqua.dtu.dk>, and Kasper Kristensen <kkr@aqua.dtu.dk>.
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//   * Neither the name of the assessment tool SAM nor the
//     names of its contributors may be used to endorse or promote products
//     derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
// ARE DISCLAIMED. IN NO EVENT SHALL ANDERS NIELSEN, CASPER BERG OR KASPER 
// KRISTENSEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  --------------------------------------------------------------------------
 
#include <TMB.hpp>
#include <iostream>

/* Calc ref point */
//using namespace TMBad;
// Utility 1: Construct slice function from objective function
template<class Type>
struct NestedOptimizer {
  const objective_function<Type>* orig;
  objective_function<TMBad::ad_aug> obj;
  TMBad::ADFun<> current_tape;
  newton::newton_config cfg;
  NestedOptimizer (const objective_function<Type>* orig,
                   newton::newton_config cfg = newton::newton_config()) :
    orig(orig), obj(*orig), cfg(cfg) {
    // Mark obj as a copy
    obj.is_copy = true;
    std::vector<TMBad::ad_aug> x(obj.theta.size()); // FIXME: DRY !!!
    for (size_t i=0; i<x.size(); i++) x[i] = obj.theta[i].Value(); // FIXME: DRY !!!

    // Create identity tape
    current_tape.glob.ad_start();

    TMBad::Independent(x);
    TMBad::Dependent(x);
    current_tape.glob.ad_stop();
  }

  std::vector<TMBad::Index> get_input_index(const char* name) {
    std::vector<TMBad::Index> input;
    for (int i=0; i<obj.thetanames.size(); i++) {
      if (!strcmp(name, obj.thetanames[i])) input.push_back(i);
    }
    return input;
  }

  std::vector<int> get_output_index(const char* name) {
    std::vector<int> output;
    SEXP names = obj.reportvector.reportnames();
    PROTECT(names);
    for (int i=0; i<LENGTH(names); i++) {
      std::cout << name << " " << CHAR(STRING_ELT(names, i)) << "\n";
      if ( ! strcmp(name, CHAR(STRING_ELT(names, i))) ) output.push_back(i);
    }
    UNPROTECT(1);
    return output;
  }
  // vector<TMBad::ad_aug> get_x(const std::vector<int> &input) {
  //   vector<TMBad::ad_aug> x(input.size());
  //   for (size_t i=0; i<input.size(); i++)
  //     x[i] = obj.theta[input[i]];
  //   return x;
  // }
  vector<TMBad::ad_aug> get_y(const std::vector<int> &output) {
    vector<TMBad::ad_aug> y(output.size());
    for (size_t i=0; i<output.size(); i++)
      y[i] = obj.reportvector.result[output[i]];
    return y;
  }
  // void set_x(const std::vector<int> &input,
  //            const vector<TMBad::ad_aug> &x) {
  //   for (size_t i=0; i<input.size(); i++)
  //     obj.theta[input[i]] = x[i];
  // }
  vector<TMBad::ad_aug> eval_obj(const vector<TMBad::ad_aug> &x,
                                 const char *name_output) {
    // Clear because we allow using same obj multiple times
    obj.index = 0;
    obj.parnames.resize(0);
    obj.reportvector.clear();
    // Eval
    for (size_t i=0; i < (size_t) x.size(); i++) obj.theta[i] = x[i];
    TMBad::ad_aug value = obj();
    if (name_output == NULL) {
      vector<TMBad::ad_aug> ans(1);
      ans[0] = value;
      return ans;
    }
    //output.resize(0);
    return get_y( get_output_index(name_output) );
  }
  // vector<TMBad::ad_aug> argmin(newton::newton_config cfg = newton::newton_config() ) {
  //   // if (output.size() != 1)
  //   //   Rf_error("Minimization requires one dimensional output");
  //   return newton::Newton(*this, get_x(), cfg);
  // }
  // void argmin_inplace(newton::newton_config cfg = newton::newton_config() ) {
  //   vector<TMBad::ad_aug> x_hat = argmin(cfg);
  //   set_x( x_hat );
  // }
  // Select output
  void set_output(const char *name) {
    if (current_tape.Range() != (size_t) obj.theta.size())
      Rf_error("Incompatible output dimension of current tape");
    TMBad::ADFun<> ans;
    //std::vector<TMBad::ad_aug> x = obj.theta;
    ans.glob.ad_start();
    
    std::vector<double> tmp = current_tape.DomainVec();
    std::vector<TMBad::ad_aug> x(tmp.begin(), tmp.end());
    // std::vector<TMBad::ad_aug> x(obj.theta.size()); // FIXME: DRY !!!
    // for (size_t i=0; i<x.size(); i++) x[i] = obj.theta[i].Value(); // FIXME: DRY !!!

    TMBad::Independent(x);
    x = current_tape(x); // Replay current tape
    std::vector<TMBad::ad_aug> y = eval_obj(x, name);
    TMBad::Dependent(y);
    ans.glob.ad_stop();
    current_tape = ans;
  }
  // argmin("F")
  void set_argmin(const char *name) {
    if (current_tape.Range() != 1)
      Rf_error("Minimization requires one dimensional output");
    TMBad::ADFun<> ans;
    //std::vector<TMBad::ad_aug> x = obj.theta;
    ans.glob.ad_start();

    std::vector<double> tmp = current_tape.DomainVec();
    std::vector<TMBad::ad_aug> x(tmp.begin(), tmp.end());
    //std::vector<TMBad::ad_aug> x(obj.theta.size()); // FIXME: DRY !!!
    //for (size_t i=0; i<x.size(); i++) x[i] = obj.theta[i].Value(); // FIXME: DRY !!!

    TMBad::Independent(x);
    //set_input(name); // FIXME: input = get_input("F")
    //std::vector<TMBad::Index> input(this->input.begin(), this->input.end());
    std::vector<TMBad::Index> input = get_input_index(name);
    newton::slice<> S(current_tape, input);
    // FIXME: WTF:
    S.x = x;
    // FIXME: newton.hpp add member 'start' to slice ???
    vector<TMBad::ad_aug> start = TMBad::subset(x, input);
    std::vector<TMBad::ad_aug> sol = newton::Newton(S, start, cfg);
    for (size_t i=0; i<input.size(); i++) x[input[i]] = sol[i];
    TMBad::Dependent(x);
    ans.glob.ad_stop();
    current_tape = ans;
  }
  void select_arg(const char *name) {
    if (current_tape.Range() != (size_t) obj.theta.size())
      Rf_error("'select_arg' requires output dimension equal to objective input dimension");
    std::vector<TMBad::Index> input = get_input_index(name);
    current_tape.glob.dep_index =
      TMBad::subset(current_tape.glob.dep_index, input);
  }
};

// NestedOptimizer objective_slice(objective_function<TMBad::ad_aug> &obj, const char *name_input, const char *name_output) {
//   NestedOptimizer s(obj);
//   s.set_input(name_input);
//   //s.set_output(name_output);
//   s.name_output = name_output;
//   return s;
// }
// Utility 2: Construct ADreport function from objective function


/* Parameter transform */
template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template <class Type> 
Type square(Type x){return x*x;}

 
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(fleetTypes); 
  DATA_VECTOR(sampleTimes);
  DATA_VECTOR(years);
  DATA_INTEGER(nobs);
  DATA_VECTOR(idx1);
  DATA_VECTOR(idx2);
  DATA_ARRAY(obs);
  DATA_ARRAY(propMat);
  DATA_ARRAY(stockMeanWeight); 
  DATA_ARRAY(catchMeanWeight);
  DATA_ARRAY(natMor);
  DATA_ARRAY(landFrac);
  DATA_ARRAY(disMeanWeight);
  DATA_ARRAY(landMeanWeight);
  DATA_ARRAY(propF);
  DATA_ARRAY(propM);
  DATA_INTEGER(minAge);
  DATA_INTEGER(maxAgePlusGroup);
  DATA_IARRAY(keyLogFsta);
  DATA_ARRAY(keyLogFpar);
  DATA_ARRAY(keyQpow);
  DATA_ARRAY(keyVarF);
  DATA_ARRAY(keyVarLogN); 
  DATA_ARRAY(keyVarObs); 
  DATA_INTEGER(stockRecruitmentModelCode);
  DATA_VECTOR(fbarRange);

  PARAMETER_VECTOR(logFpar); 
  PARAMETER_VECTOR(logQpow); 
  PARAMETER_VECTOR(logSdLogFsta); 
  PARAMETER_VECTOR(logSdLogN); 
  PARAMETER_VECTOR(logSdLogObs); 
  PARAMETER(rec_loga); 
  PARAMETER(rec_logb); 
  PARAMETER(rho); 
  PARAMETER_VECTOR(logScale); 
  PARAMETER_VECTOR(logScaleSSB); 
  PARAMETER_VECTOR(logPowSSB); 
  PARAMETER_VECTOR(logSdSSB); 
  PARAMETER_ARRAY(U);
  DATA_INTEGER(nlogF);
  DATA_INTEGER(nlogN);

  array<Type> logF(nlogF,U.cols());
  array<Type> logN(nlogN,U.cols());
  for(int i=0;i<nlogN;i++)
    for(int j=0;j<U.cols();j++)
      logN(i,j)=U(i,j);
  for(int i=0;i<nlogF;i++)
    for(int j=0;j<U.cols();j++)
      logF(i,j)=U(i+nlogN,j);


  int timeSteps=logF.dim[1];
  int stateDimF=logF.dim[0];
  int stateDimN=logN.dim[0];
  //Type rho=f(logit_rho);
  vector<Type> sdLogFsta=exp(logSdLogFsta);
  vector<Type> varLogN=exp(logSdLogN*Type(2.0));
  vector<Type> varLogObs=exp(logSdLogObs*Type(2.0));
  vector<Type> ssb(timeSteps);
  vector<Type> logssb(timeSteps);

  //First take care of F
  matrix<Type> fvar(stateDimF,stateDimF);
  matrix<Type> fcor(stateDimF,stateDimF);
  vector<Type> fsd(stateDimF);

  for(int i=0; i<stateDimF; ++i){
    for(int j=0; j<stateDimF; ++j){
      if(i!=j){fcor(i,j)=rho;}else{fcor(i,j)=1.0;}
    }
    fsd(i)=sdLogFsta(CppAD::Integer(keyVarF(0,i)));
  }
  for(int i=0; i<stateDimF; ++i){
    for(int j=0; j<stateDimF; ++j){
      fvar(i,j)=fsd(i)*fsd(j)*fcor(i,j);
    }
  }
  using namespace density;
  MVNORM_t<Type> neg_log_densityF(fvar);
  Type ans=0;
  for(int i=1;i<timeSteps;i++){    
     ans+=neg_log_densityF(logF.col(i)-logF.col(i-1)); // F-Process likelihood
     SIMULATE {
       logF.col(i) = logF.col(i-1) + neg_log_densityF.simulate();
     }
  }
 
  for(int i=0;i<timeSteps;i++){ // calc ssb
    ssb(i)=0.0;    
    for(int j=0; j<stateDimN; ++j){
      ssb(i)+=exp(logN(j,i))*exp(-exp(logF((keyLogFsta(0,j)),i))*propF(i,j)-natMor(i,j)*propM(i,j))*propMat(i,j)*stockMeanWeight(i,j);
    }
    logssb(i)=log(ssb(i));
  }

  //Now take care of N
  matrix<Type> nvar(stateDimN,stateDimN);
  for(int i=0; i<stateDimN; ++i){
    for(int j=0; j<stateDimN; ++j){
      if(i!=j){nvar(i,j)=0.0;}else{nvar(i,j)=varLogN(CppAD::Integer(keyVarLogN(0,i)));}
    }
  }
  MVNORM_t<Type> neg_log_densityN(nvar);
  vector<Type> predN(stateDimN); 
  for(int i=1;i<timeSteps;i++){ 
    if(stockRecruitmentModelCode==0){ // straight RW 
      predN(0)=logN(0,i-1);
    }else{
      if(stockRecruitmentModelCode==1){//ricker
        predN(0)=rec_loga+log(ssb(i-1))-exp(rec_logb)*ssb(i-1); 
      }else{
        if(stockRecruitmentModelCode==2){//BH
          predN(0)=rec_loga+log(ssb(i-1))-log(1.0+exp(rec_logb)*ssb(i-1)); 
        }else{
          error("SR model code not recognized");
        }
      }
    }

    for(int j=1; j<stateDimN; ++j){
      predN(j)=logN(j-1,i-1)-exp(logF((keyLogFsta(0,j-1)),i-1))-natMor(i-1,j-1); 
    }  
    if(maxAgePlusGroup==1){
      predN(stateDimN-1)=log(exp(logN(stateDimN-2,i-1)-exp(logF((keyLogFsta(0,stateDimN-2)),i-1))-natMor(i-1,stateDimN-2))+
                             exp(logN(stateDimN-1,i-1)-exp(logF((keyLogFsta(0,stateDimN-1)),i-1))-natMor(i-1,stateDimN-1))); 
    }
    ans+=neg_log_densityN(logN.col(i)-predN); // N-Process likelihood
    SIMULATE {
      logN.col(i) = predN + neg_log_densityN.simulate();
    }
  }


  // Now finally match to observations
  int f, ft, a, y; 
  int minYear=CppAD::Integer((obs(0,0)));
  Type predObs=0, zz, var;
  for(int i=0;i<nobs;i++){
    y=CppAD::Integer(obs(i,0))-minYear;
    f=CppAD::Integer(obs(i,1));
    ft=CppAD::Integer(fleetTypes(f-1));
    a=CppAD::Integer(obs(i,2))-minAge;
    zz=exp(logF((keyLogFsta(0,a)),y))+natMor(y,a);
    
    if(ft==0){// residual fleet
      predObs=logN(a,y)-log(zz)+log(1-exp(-zz));
      if((keyLogFsta(f-1,a))>(-1)){
        predObs+=logF((keyLogFsta(0,a)),y);
      }
    }else{
      if(ft==1){// comm fleet
        //cerr<<"Not implemented yet!!!"<<endl;  
      }else{
        if(ft==2){// survey
          predObs=logN(a,y)-zz*sampleTimes(f-1);
          if(CppAD::Integer(keyQpow(f-1,a))>(-1)){
            predObs*=exp(logQpow(CppAD::Integer(keyQpow(f-1,a)))); 
          }
          if(CppAD::Integer(keyLogFpar(f-1,a))>(-1)){
            predObs+=logFpar(CppAD::Integer(keyLogFpar(f-1,a)));
          }
        }else{
          if(ft==3){// SSB survey -- nevermind for now 
            //cerr<<"Not implemented yet!!!"<<endl;  
          }else{
            if(ft==4){// SSB survey -- nevermind for now 
              //cerr<<"Not implemented yet!!!"<<endl;  
            }
          }
        } 
      }
    }      
    var=varLogObs(CppAD::Integer(keyVarObs(f-1,a)));
    ans+=-dnorm(log(obs(i,3)),predObs,sqrt(var),true);
    SIMULATE {
      obs(i,3) = exp( rnorm(predObs, sqrt(var)) ) ;
    }
  }

  SIMULATE {
    REPORT(logF);
    REPORT(logN);
    REPORT(obs);
  }

  // Some objective
  Type SSB = rho*rho + U[0]*U[0];
  ADREPORT(ans);
  ADREPORT(SSB);
  // Avoid infinite recursion
  DATA_INTEGER(flag);
  if ( flag && ! this -> is_copy ) {
    newton::newton_config cfg = newton::newton_config();
    cfg.trace = true;
    cfg.sparse = true;
    NestedOptimizer<Type> Nopt(this, cfg);
    Nopt.set_output("ans");
    Nopt.set_argmin("U");
    Nopt.set_output("SSB");
    Nopt.set_argmin("rho");
    Nopt.select_arg("rho");
    // Replace U by Uhat in 'obj'
    // objective_slice(obj, "U", "ans").argmin_inplace();
    // vector<TMBad::ad_aug> Fhat = objective_slice(obj, "F", "SSB").argmin();
    // REPORT(Fhat);
    vector<Type> rho_hat = Nopt.current_tape(this->theta);
    REPORT(rho_hat);
    ADREPORT(rho_hat);
  }
  // ADREPORT(logN);
  // ADREPORT(logF);
  return ans;
}
