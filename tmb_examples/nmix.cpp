// nmix example from https://groups.nceas.ucsb.edu/non-linear-modeling/projects/nmix
#include <TMB.hpp>

template<class Type>
Type nll_group(int i, Type p0,Type p1,Type log_lambda, Type log_sigma, 
	       vector<Type> u,matrix<Type> y,vector<Type> N, vector<Type> x, 
	       matrix<Type> IDind){
  using CppAD::Integer;
  int T=y.cols();
  int S=N.size();
  Type sigma = exp(log_sigma);
  Type lambda = exp(log_lambda);
  Type e=1e-12;
  Type nll=0;
  vector<Type> logf(S);
  vector<Type> logg(S);
  vector<Type> fg(S);
  vector<Type> tmp=Type(-1.0)*(p0 + p1*x(i) + sigma*u);
  vector<Type> p = Type(1.0)/(Type(1.0)+exp(tmp));
  for(int k=0;k<S;k++){
    logf(k) = log_lambda*N(k) - (lambda + lgamma(N(k)+1));
  }
  Type tmp1,tmp2,tmp3;
  for(int k=0;k<S;k++) {
    logg(k) = 0;
    for(int j=0;j<T;j++) {
      if(N(k)>=y(i,j)){
	tmp1=lgamma(N(k)+1) - lgamma(y(i,j)+1) - lgamma(N(k)-y(i,j)+1);
	tmp2=log(p(Integer(IDind(i,j)))+e)*y(i,j);
	tmp3=log(Type(1.0) + e - p(Integer(IDind(i,j))))*( N(k)-y(i,j));
	logg(k)+=tmp1+tmp2+tmp3;
      }
      else {
	logg(k) = -1000;
      }
    }
  }
  fg=exp(logf+logg);
  nll -= log(e + sum(fg));
  return nll;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* data section */
  DATA_INTEGER(R);              // Number of sites
  DATA_VECTOR(N);               // Possible values of N
  DATA_VECTOR(nID);             // Number of observers present at a site
  DATA_FACTOR(ID);              // IDs of observer present at each site
  DATA_MATRIX(IDind);           // Group ID RT matrix
  DATA_FACTOR(IDfac);           // To split ID (instead of ragged matrix)
  DATA_VECTOR(x);               // Site-specific covariate
  DATA_MATRIX(y);               // Count data with R rows and T columns
  /* Parameter section */
  PARAMETER(log_lambda);
  PARAMETER(p0);
  PARAMETER(p1);
  PARAMETER(log_sigma);         // log of random effect SD
  PARAMETER_VECTOR(u);          // Length nG
  vector<vector<int> > idspl=split(ID,IDfac);
  /* Procedure section */
  Type nll=0;
  nll+=Type(.5)*(u*u).sum();
  for(int i=0;i<R;i++){
    nll+=nll_group(i, p0,p1,log_lambda,log_sigma,u(idspl(i)),y,N,x,IDind);
  }
  return nll;
}
