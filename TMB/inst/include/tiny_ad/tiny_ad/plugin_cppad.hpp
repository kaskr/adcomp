/* The following methods are required for a CppAD Base type, see the
   file 'cppad/local/base_double.hpp'.

*/
template<class T, class V>
ad<T, V> abs (const ad<T, V> &x){
  return fabs(x);
}
template<class T>
T CondExpOp( 
            enum CppAD::CompareOp cop          ,
            const T&       left         ,
            const T&       right        , 
            const T&       exp_if_true  , 
            const T&       exp_if_false )
{	return CppAD::CondExpTemplate(cop, left, right, exp_if_true, exp_if_false);
}
template<class T>  int Integer(const T& x)
{	return Integer(x.value); }
template<>  int Integer(const double& x)
{	return static_cast<int>(x); }
template<class T>  bool GreaterThanZero(const T& x)
{	return x > 0.; }
template<class T>  bool GreaterThanOrZero(const T& x)
{	return x >= 0.; }
template<class T>  bool LessThanZero(const T& x)
{	return x < 0.; }
template<class T>  bool LessThanOrZero(const T& x)
{	return x <= 0.; }
template<class T>  bool abs_geq(const T& x, const T& y)
{	return fabs(x) >= fabs(y); }
template<class T> bool IdenticalPar(const T& x)
{	return true; }
template<class T> bool IdenticalZero(const T& x)
{	return (x == 0.); }
template<class T> bool IdenticalOne(const T& x)
{	return (x == 1.); }
template<class T> bool IdenticalEqualPar(const T& x, const T& y)
{	return (x == y); }
