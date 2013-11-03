/** 
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
   \brief Namespace to construct multivariate Gaussian distributions via C++ templates

   A particular multivariate normal distribution is implemented as a templated C++ class. 
   Let us take the generic zero-mean multivariate normal distribution <tt>MVNORM_t</tt> with covariance matrix 
   <tt>Sigma</tt> as an example. The <tt>_t</tt> symbol attached to the class name reminds us that we are dealing with a class
   (of a particular C++ type). There are two operations that we can do on objects from the <tt>MVNORM_t</tt> class:
   <list>
   <item> Declare and initialize in terms of one or more parameters (e.g. <tt>Sigma</tt>)     </item>
   <item> Evaluate the negative log-likelihood density at specified point (e.g. a vector <tt>u</tt>)   </item>   
   </list>
   An example is
   \code
     // Sigma and u are objects that have allready been defined
     MVNORM_t<Type> neg_log_density(Sigma);	// Create object from covariance matrix Sigma
     Type ans = neg_log_density(u);		// Evaluate neg. log-likelihood
   \endcode
   Comments:
   <list> 
   <item> The template argument <tt><Type></tt> must always be included (as in many other places in TMB) but can be ignored
	   from a user perspective.</item>
   <item> The object vi define here is called <tt>neg_log_density</tt>. You can choose whatever name you like,
	  but <tt>neg_log_density</tt> reminds you that the only thing you will do with it is to 
	  evaluate the negative log-likelihood. </item>
   <item> The dimensions of <tt>Sigma</tt> and <tt>u</tt> must match. </item>   
   </list>

   New classes of distributions can be built recursively from existing distributions. The behaviour of the different
   distributions differ in this respect. A very useful example is the
   Kronecker product (http://en.wikipedia.org/wiki/Kronecker_product) of two multivariate normal distributions (see detailed description of <tt>SEPARABLE_t</tt>):
    \code
     MVNORM_t<Type> Gauss1(Sigma1);	// First normal distribution (with covariance Sigma1)
     MVNORM_t<Type> Gauss2(Sigma2);	// Second normal distribution (with covariance Sigma2)
     SEPARABLE_t<MVNORM_t<Type> , MVNORM_t<Type> > Gauss3(Gauss1,Gauss2);
     SEPARABLE_t<MVNORM_t<Type>, MVNORM_t<Type>> Gauss3(Gauss1,Gauss2);
     Type ans = Gauss3(u);		// Evaluate neg. log-likelihood
    \endcode
    where <tt>u</tt> must be of appropriate dimension.
*/
namespace density{
  using namespace tmbutils;
#include "density.cpp"
} // End namespace
