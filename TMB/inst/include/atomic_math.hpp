namespace atomic {
/* Namespace with double versions of R special math library */
namespace Rmath {
  #include <Rmath.h>
  // Macros do not respect the namespace limits.
  // Some of them will conflict with other TMB functions.
  #undef dnorm
  #undef pnorm
  #undef qnorm

  #include <R_ext/Applic.h>
  void integrand_D_incpl_gamma_shape(double *x, int n, void *ex){
    double* parms=(double*)ex;
    for(int i=0;i<n;i++) x[i] = exp(-x[i]) * pow(x[i],parms[0]-1.0) * pow(log(x[i]),parms[1]);
  }
  /* n'th order derivative of incomplete gamma wrt. shape parameter */
  double D_incpl_gamma_shape(double x, double shape, double n){
    double a=0;
    double b=x;
    double epsabs=1e-8;
    double epsrel=1e-8;
    double result=0;
    double abserr=10000;
    int neval=10000;
    int ier=0;
    int limit=100;
    int lenw = 4 * limit;
    int last=0;
    int* iwork =  Calloc(limit, int);
    double* work = Calloc(lenw, double);
    double ex[2];
    ex[0]=shape;
    ex[1]=n;
    Rdqags(integrand_D_incpl_gamma_shape, ex, &a, &b,
	   &epsabs, &epsrel,
	   &result, &abserr, &neval, &ier,
	   &limit, &lenw, &last, iwork, work);
    Free(iwork);
    Free(work);
    if(ier!=0)warning("Integrate incomplete gamma function unreliable");
    return result;
  }

  double inv_incpl_gamma(double y, double shape){
    double p=y/gammafn(shape);
    double scale=1.0;
    return qgamma(p, shape, scale, 1, 0);
  }

  /* n'th order derivative of log gamma function */
  double D_lgamma(double x, double n){
    if(n<.5)return lgammafn(x);
    else return psigamma(x,n-1.0);
  }

}

#include "atomic_macro.hpp"

/* Utilities for conversion between CppAD::vector and Eigen::Matrix */
template<class Type>
matrix<Type> vec2mat(CppAD::vector<Type> x, int m, int n, int offset=0){
  matrix<Type> res(m,n);
  for(int i=0;i<m*n;i++)res(i)=x[i+offset];
  return res;
}
template<class Type>
CppAD::vector<Type> mat2vec(matrix<Type> x){
  int n=x.size();
  CppAD::vector<Type> res(n);
  for(int i=0;i<n;i++)res[i]=x(i);
  return res;
}

/* Derivative of pnorm1 (based on functions with known derivatives) */
template<class Type>
Type dnorm1(Type x){
  return Type(1.0/sqrt(2.0*M_PI)) * exp(-Type(.5)*x*x);
}

TMB_ATOMIC_VECTOR_FUNCTION(
			   // ATOMIC_NAME
			   pnorm1
			   ,
			   // OUTPUT_DIM
			   1
			   ,
			   // ATOMIC_DOUBLE
			   ty[0] = Rmath::pnorm5(tx[0],0,1,1,0);
			   ,
			   // ATOMIC_REVERSE
			   px[0] = dnorm1(tx[0]) * py[0];
			   )

TMB_ATOMIC_VECTOR_FUNCTION(
			   // ATOMIC_NAME
			   qnorm1
			   ,
			   // OUTPUT_DIM
			   1,
			   // ATOMIC_DOUBLE
			   ty[0] = Rmath::qnorm5(tx[0],0,1,1,0);
			   ,
			   // ATOMIC_REVERSE
			   px[0] = Type(1) / dnorm1(ty[0]) * py[0];
			   )

TMB_ATOMIC_VECTOR_FUNCTION(
			   // ATOMIC_NAME
			   D_incpl_gamma_shape
			   ,
			   // OUTPUT_DIM
			   1
			   ,
			   // ATOMIC_DOUBLE
			   ty[0]=Rmath::D_incpl_gamma_shape(tx[0],tx[1],tx[2]);
			   ,
			   // ATOMIC_REVERSE
			   px[0] = exp(-tx[0])*pow(tx[0],tx[1]-Type(1.0))*pow(log(tx[0]),tx[2]) * py[0];
			   CppAD::vector<Type> tx_(tx);
			   tx_[2] = tx_[2] + Type(1.0);  // Add one to get partial wrt. tx[1]
			   px[1] = D_incpl_gamma_shape(tx_)[0] * py[0];
			   px[2] = Type(0);
			   )

TMB_ATOMIC_VECTOR_FUNCTION(
			   // ATOMIC_NAME
			   inv_incpl_gamma
			   ,
			   // OUTPUT_DIM
			   1
			   ,
			   // ATOMIC_DOUBLE
			   ty[0]=Rmath::inv_incpl_gamma(tx[0],tx[1]);
			   ,
			   // ATOMIC_REVERSE
			   Type value = ty[0];
			   Type shape = tx[1];
			   Type tmp = exp(-value)*pow(value,shape-Type(1));
			   px[0] = 1.0 / tmp * py[0];
			   CppAD::vector<Type> arg(3);
			   arg[0] = value;
			   arg[1] = shape;
			   arg[2] = Type(1); // 1st order partial
			   px[1] = -D_incpl_gamma_shape(arg)[0] / tmp * py[0];
			   )

TMB_ATOMIC_VECTOR_FUNCTION(
			   // ATOMIC_NAME
			   D_lgamma
			   ,
			   // OUTPUT_DIM
			   1
			   ,
			   // ATOMIC_DOUBLE
			   ty[0]=Rmath::D_lgamma(tx[0],tx[1]);
			   ,
			   // ATOMIC_REVERSE
			   CppAD::vector<Type> tx_(2);
			   tx_[0]=tx[0];
			   tx_[1]=tx[1]+Type(1.0);
			   px[0] = D_lgamma(tx_)[0] * py[0];
			   px[1] = Type(0);
			   )

template<class Type> /* Header of matmul interface */
matrix<Type> matmul(matrix<Type> x, matrix<Type> y);
TMB_ATOMIC_VECTOR_FUNCTION(
			   // ATOMIC_NAME
			   matmul
			   ,
			   // OUTPUT_DIM (currently square matrix only)
			   tx.size()/2
			   ,
			   // ATOMIC_DOUBLE
			   int n=sqrt(tx.size()/2);
			   matrix<double> X = vec2mat(tx, n, n, 0);
			   matrix<double> Y = vec2mat(tx, n, n, n*n);
			   matrix<double> res = X * Y;       // Use Eigen matrix multiply
			   for(int i=0;i<n*n;i++)ty[i] = res(i);
			   ,
			   // ATOMIC_REVERSE (W*Y^T, X^T*W)
			   int n=sqrt(ty.size());
			   matrix<Type> Xt = vec2mat(tx,n,n,0).transpose();
			   matrix<Type> Yt = vec2mat(tx,n,n,n*n).transpose();
			   matrix<Type> W = vec2mat(py,n,n);
			   matrix<Type> res1 = matmul(W, Yt); // W*Y^T
			   matrix<Type> res2 = matmul(Xt, W); // X^T*W
			   for(int i=0;i<n*n;i++){
			     px[i]     = res1(i);
			     px[i+n*n] = res2(i);
			   }
			   )

TMB_ATOMIC_VECTOR_FUNCTION(
			   // ATOMIC_NAME
			   matinv
			   ,
			   // OUTPUT_DIM
			   tx.size()
			   ,
			   // ATOMIC_DOUBLE
			   int n=sqrt(tx.size());
			   matrix<double> X(n,n);
			   for(int i=0;i<n*n;i++){X(i)=tx[i];}
			   matrix<double> res=X.inverse();   // Use Eigen matrix inverse (LU)
			   for(int i=0;i<n*n;i++)ty[i]=res(i);
			   ,
			   // ATOMIC_REVERSE  (-f(X)^T*W*f(X)^T)
			   int n=sqrt(ty.size());
			   matrix<Type> W=vec2mat(py,n,n);   // Range direction
			   matrix<Type> Y=vec2mat(ty,n,n);   // f(X)
			   matrix<Type> Yt=Y.transpose();    // f(X)^T
			   matrix<Type> tmp=matmul(W,Yt);    // W*f(X)^T
			   matrix<Type> res=-matmul(Yt,tmp); // -f(X)^T*W*f(X)^T
			   px=mat2vec(res);
			   )

TMB_ATOMIC_VECTOR_FUNCTION(
			   // ATOMIC_NAME
			   logdet
			   ,
			   // OUTPUT_DIM
			   1
			   ,
			   // ATOMIC_DOUBLE
			   int n=sqrt(tx.size());
			   matrix<double> X=vec2mat(tx,n,n);
			   matrix<double> LU=X.lu().matrixLU();    // Use Eigen LU decomposition
			   vector<double> LUdiag = LU.diagonal();
			   double res=LUdiag.abs().log().sum();    // TODO: currently PD only - take care of sign.
			   ty[0] = res;
			   ,
			   // ATOMIC_REVERSE  (X^-1*W[0])
			   CppAD::vector<Type> invX=matinv(tx);
			   for(int i=0;i<tx.size();i++)px[i]=invX[i]*py[0];
			   )

/* ================================== INTERFACES
*/

template<class Type>
matrix<Type> matmul(matrix<Type> x, matrix<Type> y){
  CppAD::vector<Type> arg(x.size()+y.size());
  for(int i=0;i<x.size();i++){arg[i]=x(i);}
  for(int i=0;i<y.size();i++){arg[i+x.size()]=y(i);}
  return vec2mat(matmul(arg),x.rows(),y.cols());
}

template<class Type>
matrix<Type> matinv(matrix<Type> x){
  int n=x.rows();
  return vec2mat(matinv(mat2vec(x)),n,n);
}

template<class Type>
Type logdet(matrix<Type> x){
  return logdet(mat2vec(x))[0];
}

template<class Type>
Type pnorm(Type q, Type mean = 0, Type sd = 1){
  CppAD::vector<Type> tx(1);
  tx[0] = (q - mean) / sd;
  return pnorm1(tx)[0];
}

template<class Type>
Type qnorm(Type p, Type mean = 0, Type sd = 1){
  CppAD::vector<Type> tx(1);
  tx[0] = p;
  return sd*qnorm1(tx)[0] + mean;
}

template<class Type>
Type lgamma(Type x){
  CppAD::vector<Type> tx(2);
  tx[0] = x;
  tx[1] = Type(0);
  return D_lgamma(tx)[0];
}

template<class Type>
Type pgamma(Type q, Type shape, Type scale = 1){
  CppAD::vector<Type> tx(3);
  tx[0] = q/scale;
  tx[1] = shape;
  tx[2] = Type(0); // 0'order deriv
  return D_incpl_gamma_shape(tx)[0] / exp(lgamma(shape));
}

template<class Type>
Type qgamma(Type q, Type shape, Type scale = 1){
  CppAD::vector<Type> tx(2);
  tx[0] = q * exp(lgamma(shape));
  tx[1] = shape;
  return inv_incpl_gamma(tx)[0] * scale;
}

/* Temporary test of dmvnorm implementation based on atomic symbols.
   Should reduce tape size from O(n^3) to O(n^2).
*/
template<class Type>
Type nldmvnorm(vector<Type> x, matrix<Type> Sigma){
  matrix<Type> Q=matinv(Sigma);
  Type logdetQ = -logdet(Sigma);
  Type quadform = (x*(Q*x)).sum();
  return -Type(.5)*logdetQ + Type(.5)*quadform + x.size()*Type(log(sqrt(2.0*M_PI)));
}

} /* End namespace atomic */

/* User defined atomic functions */
#include "checkpoint_macro.hpp"
