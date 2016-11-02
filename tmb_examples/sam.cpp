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
  // ADREPORT(logN);
  // ADREPORT(logF);
  return ans;
}
