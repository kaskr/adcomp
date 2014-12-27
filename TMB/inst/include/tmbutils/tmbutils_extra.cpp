namespace tmbutils{
// Convenience utilites
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

} // End namespace

/** 
   \brief Collection of multivariate Gaussian distributions (members listed in \ref density.cpp)

   \ingroup Densities
   
   For use of the namespace see \ref Densities   
*/
namespace density{
  using namespace tmbutils;
#include "density.cpp"
} // End namespace
