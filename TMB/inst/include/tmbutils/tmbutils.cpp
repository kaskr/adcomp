/** \file 
   \brief Namespace of utility functions for TMB
*/
namespace tmbutils{
//#include <Eigen/Dense>
using namespace Eigen;
#include "vector.cpp"
#include "array.cpp"
#include "spmat.cpp"

#include "kronecker.cpp"
#include "matexp.cpp"
#include "splines.cpp"
#include "order.cpp"

template<class Type, class T1, class T2>
vector<Type> dnorm(vector<Type> x, T1 mean, T2 sd, int give_log=0)
{
  vector<Type> logres;
  x=(x-mean)/sd;
  logres=-log(Type(sqrt(2*M_PI))*sd)-Type(.5)*x*x;
  if(give_log)return logres; else return exp(logres);
}
  
template <class Type, class From>
vector<Type> asVector(From *px, int n){
  vector<Type> x(n);
  for(int i=0;i<n;i++)x[i]=Type(px[i]);
  return x;
}

#if defined(R_R_H) 
template <class Type>
array<Type> asArray(SEXP x)
{
  if(!isArray(x))error("NOT AN ARRAY!");
  SEXP dim=getAttrib(x,R_DimSymbol);
  vector<int> d=asVector<int,int>(INTEGER(dim), LENGTH(dim));
  vector<Type> y=asVector<Type,double>(REAL(x), LENGTH(x));
  return array<Type>(y,d);
}


#endif


} // End namespace

/** 
   \brief Collection of multivariate Gaussian distributions (members listed in \ref density.cpp)

   \ingroup Densities

   The purpose of the name space is to provide multivariate normal distributions
   useful for classical multivariate analysis, time series, spatial model, and space-time models.
   Let us take the generic zero-mean multivariate normal distribution <tt>MVNORM</tt> 
   with covariance matrix \c Sigma as an example. To evaluate the negative log-likelihood
   at the point <tt>u</tt> you simply write
   \code
     PARAMETER_VECTOR(u);				// Random vector
     PARAMETER_MATRIX(Sigma)			// Covariance matrix
     Type ans = MVNORM(Sigma)(u);		// Evaluate neg. log-likelihood
   \endcode
   
   You will also see <tt>MVNORM_t</tt> (note the <tt>_t</tt> extention) which is a C++ class 
   that implements the distribution. You will use <tt>MVNORM_t</tt> to build more complicated
   distributions, for instance via the Kronecker product (http://en.wikipedia.org/wiki/Kronecker_product).
   Now, let us define two distributions:
	
    \code
     MVNORM_t<Type> Gauss1(Sigma1);	     // Define a MVNORM called "Gauss1"
     MVNORM_t<Type> Gauss2(Sigma2);	     // Define a MVNORM called "Gauss2"
    \endcode
  i.e. we have got two multivariate normal distributions called \c Gauss1 and \c Gauss2. If we want
  to evaluate these we can use
   \code
     PARAMETER_VECTOR(u);
     Type ans = Gauss1(Sigma)(u);		// Evaluate neg. log-likelihood
     ans += Gauss2(Sigma)(u);		// Evaluate neg. log-likelihood
   \endcode
  but that is not our purpose here. We want to create a new distribution <tt>Gauss3</tt>
  which is the Kronecker product of <tt>Gauss1</tt> and <tt>Gauss2</tt>
    \code  
     PARAMETER_VECTOR(v);		// dim(v) = dim(u)^2
     SEPARABLE_t<MVNORM_t<Type>, MVNORM_t<Type>> Gauss3(Gauss1,Gauss2);
     Type ans = Gauss3(u2);		// Evaluate neg. log-likelihood
    \endcode
    This is complicated stuff, but shows the power of TMB!

    To be precise, \c Gauss 3 is a zero-mean multivariate normal distribution with
    covariance matrix \c Sigma3 given by the pseudo-code
    \code
     Sigma3 = inv(kronecker(inv(Sigma1),inv(Sigma2)))
    \endcode
       
   For more details about Kronecker products in TMB see <tt>SEPARABLE_t</tt>.
*/
namespace density{
  using namespace tmbutils;
#include "density.cpp"
} // End namespace
