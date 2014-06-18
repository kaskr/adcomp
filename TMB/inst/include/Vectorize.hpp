
/* Macros to do vectorization */
#define VECTORIZE1(FUN)				\
template <class Type>				\
vector<Type> FUN(const vector<Type> &x)		\
{						\
  vector<Type> res(x.size());			\
  for(int i=0;i<x.size();i++)res[i]=FUN(x[i]);	\
  return res;					\
}						\
template <class Type>				\
matrix<Type> FUN(const matrix<Type> &x)		\
{						\
  matrix<Type> res(x.rows(),x.cols());	\
  for(int i=0;i<x.rows();i++)  		\
    for(int j=0;j<x.cols();j++)		\
      res(i,j)=FUN(x(i,j));			\
  return res;					\
}
#define VECTORIZE32(FUN)						\
template <class T,class T4,class T5>					\
vector<T> FUN(const vector<T> &x1, const vector<T> &x2, 		\
	      const vector<T> &x3, T4 x4, T5 x5){			\
  int n1=x1.size();int n2=x2.size();int n3=x3.size();			\
  int n=fmax(fmax(n1,n2),n3);						\
  vector<T> res(n);							\
  for(int i=0;i<n;i++)res[i]=FUN(x1[i%n1],x2[i%n2],x3[i%n3],x4,x5);	\
  return res;								\
}
#define VECTORIZE31(FUN)						\
template <class T>							\
vector<T> FUN(const vector<T> &x1, const vector<T> &x2,			\
	      const vector<T> &x3, int x4){				\
  int n1=x1.size();int n2=x2.size();int n3=x3.size();			\
  int n=int(fmax(fmax(n1,n2),n3));					\
  vector<T> res(n);							\
  for(int i=0;i<n;i++)res[i]=FUN(x1[i%n1],x2[i%n2],x3[i%n3],x4);	\
  return res;								\
}
#define VECTORIZE21(FUN)						\
template <class T>							\
vector<T> FUN(const vector<T> &x1, const vector<T> &x2, int x4){	\
  int n1=x1.size();int n2=x2.size();					\
  int n=int(fmax(n1,n2));						\
  vector<T> res(n);							\
  for(int i=0;i<n;i++)res[i]=FUN(x1[i%n1],x2[i%n2],x4);			\
  return res;								\
}

#define VECTORIZE2(FUN)					\
template <class T>					\
vector<T> FUN(const vector<T> &x1, const T &x2){	\
  int n=x1.size();					\
  vector<T> res(n);					\
  for(int i=0;i<n;i++)res[i]=FUN(x1[i],x2);		\
  return res;						\
}

/**	\brief Vectorization macro 3.Tti
	
	For three-arguments functions (Type, Type, int). Vectorize first argument.
	*/
#define VECTORIZE3_Tti(FUN)							\
template <class Type>								\
vector<Type> FUN(const vector<Type> &arg1, Type arg2, int arg3)			\
{										\
	vector<Type> res(arg1.size());						\
	for(int i=0;i<arg1.size();i++) res[i] = pexp(arg1[i],arg2,arg3);	\
	return res;								\
}

/**	\brief Vectorization macro 3.tTi
	
	For three-arguments functions (Type, Type, int). Vectorize second argument.
	*/
#define VECTORIZE3_tTi(FUN)							\
template <class Type>								\
vector<Type> FUN(Type arg1, const vector<Type> &arg2, int arg3)			\
{										\
	vector<Type> res(arg2.size());						\
	for(int i=0;i<arg2.size();i++) res[i] = pexp(arg1,arg2[i],arg3);	\
	return res;								\
}

/**	\brief Vectorization macro 3.TTi
	
	For three-arguments functions (Type, Type, int). Vectorize first and second arguments.
	*/
#define VECTORIZE3_TTi(FUN)								\
template <class Type>									\
vector<Type> FUN(const vector<Type> &arg1, const vector<Type> &arg2, int arg3)		\
{											\
	int n1 = arg1.size();								\
	int n2 = arg2.size();								\
	int n = fmax(n1,n2);								\
	vector<Type> res(n);								\
	for(int i=0;i<n;i++) res[i] = pexp(arg1[i],arg2[i],arg3);			\
	return res;									\
}

/**	\brief Vectorization macro 4.Iiti
	
	For four-arguments functions (int, int, Type, int). Vectorize first argument.
	*/
#define VECTORIZE4_Iiti(FUN)							\
template <class Type>								\
vector<Type> FUN(const vector<int> &arg1, int arg2, Type arg3, int arg4)	\
{										\
	  vector<Type> res(arg1.size());					\
	  for(int i=0;i<arg1.size();i++) res[i]=FUN(arg1[i],arg2,arg3,arg4);	\
	  return res;								\
}

/**	\brief Vectorization macro 4.iIti
	
	For four-arguments functions (int, int, Type, int). Vectorize second argument.
	*/
#define VECTORIZE4_iIti(FUN)							\
template <class Type>								\
vector<Type> FUN(int arg1, const vector<int> &arg2, Type arg3, int arg4)	\
{										\
	vector<Type> res(arg2.size());						\
	for(int i=0;i<arg2.size();i++) res[i] = FUN(arg1,arg2[i],arg3,arg4);	\
	return res;								\
}											

/**	\brief Vectorization macro 4.iiTi
	
	For four-arguments functions (int, int, Type, int). Vectorize third argument.
	*/
#define VECTORIZE4_iiTi(FUN)							\
template <class Type>								\
vector<Type> FUN(int arg1, int arg2, const vector<Type> &arg3, int arg4)	\
{										\
	vector<Type> res(arg3.size());						\
	for(int i=0;i<arg3.size();i++) res[i] = FUN(arg1,arg2,arg3[i],arg4);	\
	return res;								\
}

/**	\brief Vectorization macro 4.iITi
	
	For four-arguments functions (int, int, Type, int). Vectorize second and third arguments.
	*/
#define VECTORIZE4_iITi(FUN)								\
template <class Type>									\
vector<Type> FUN(int arg1, const vector<int> &arg2, const vector<Type> &arg3, int arg4)	\
{											\
	int n2 = arg2.size();								\
	int n3 = arg3.size();								\
	int n = fmax(n2,n3);								\
	vector<Type> res(n);								\
	for(int i=0;i<n;i++) res[i] = FUN(arg1,arg2[i],arg3[i],arg4);			\
	return res;									\
}

/**	\brief Vectorization macro 4.IiTi
	
	For four-arguments functions (int, int, Type, int). Vectorize first and third arguments.
	*/
#define VECTORIZE4_IiTi(FUN)								\
template <class Type>									\
vector<Type> FUN(const vector<int> &arg1, int arg2, const vector<Type> &arg3, int arg4)	\
{											\
	int n1 = arg1.size();								\
	int n3 = arg3.size();								\
	int n = fmax(n1,n3);								\
	vector<Type> res(n);								\
	for(int i=0;i<n;i++) res[i] = FUN(arg1[i],arg2,arg3[i],arg4);			\
	return res;									\
}

/**	\brief Vectorization macro 4.IIti
	
	For four-arguments functions (int, int, Type, int). Vectorize first and second arguments.
	*/
#define VECTORIZE4_IIti(FUN)								\
template <class Type>									\
vector<Type> FUN(const vector<int> &arg1, const vector<int> &arg2, Type arg3, int arg4)	\
{											\
	int n1 = arg1.size();								\
	int n2 = arg2.size();								\
	int n = fmax(n1,n2);								\
	vector<Type> res(n);								\
	for(int i=0;i<n;i++) res[i] = FUN(arg1[i],arg2[i],arg3,arg4);			\
	return res;									\
}

/**	\brief Vectorization macro 4.IITi
	
	For four-arguments functions (int, int, Type, int). Vectorize first, second and third arguments.
	*/
#define VECTORIZE4_IITi(FUN)										\
template <class Type>											\
vector<Type> FUN(const vector<int> &arg1, const vector<int> &arg2, const vector<Type> &arg3, int arg4)	\
{													\
	int n1 = arg1.size();										\
	int n2 = arg2.size();										\
	int n3 = arg3.size();										\
	int n = fmax(fmax(n1,n2),n3);									\
	vector<Type> res(n);										\
	for(int i=0;i<n;i++) res[i] = FUN(arg1[i],arg2[i],arg3[i],arg4);				\
	return res;											\
}

/**	\brief Vectorization macro 4.Ttti
	
	For four-arguments functions (Type, Type, Type, int). Vectorize first argument.
	*/
#define VECTORIZE4_Ttti(FUN)							\
template <class Type>								\
vector<Type> FUN(const vector<Type> &arg1, Type arg2, Type arg3, int arg4)	\
{										\
	  vector<Type> res(arg1.size());					\
	  for(int i=0;i<arg1.size();i++) res[i]=FUN(arg1[i],arg2,arg3,arg4);	\
	  return res;								\
}

/**	\brief Vectorization macro 4.tTti
	
	For four-arguments functions (Type, Type, Type, int). Vectorize second argument.
	*/
#define VECTORIZE4_tTti(FUN)							\
template <class Type>								\
vector<Type> FUN(Type arg1, const vector<Type> &arg2, Type arg3, int arg4)	\
{										\
	vector<Type> res(arg2.size());						\
	for(int i=0;i<arg2.size();i++) res[i] = FUN(arg1,arg2[i],arg3,arg4);	\
	return res;								\
}											

/**	\brief Vectorization macro 4.ttTi
	
	For four-arguments functions (Type, Type, Type, int). Vectorize third argument.
	*/						\
#define VECTORIZE4_ttTi(FUN)							\
template <class Type>								\
vector<Type> FUN(Type arg1, Type arg2, const vector<Type> &arg3, int arg4)	\
{										\
	vector<Type> res(arg3.size());						\
	for(int i=0;i<arg3.size();i++) res[i] = FUN(arg1,arg2,arg3[i],arg4);	\
	return res;								\
}

/**	\brief Vectorization macro 4.tTTi
	
	For four-arguments functions (Type, Type, Type, int). Vectorize second and third arguments.
	*/
#define VECTORIZE4_tTTi(FUN)									\
template <class Type>										\
vector<Type> FUN(Type arg1, const vector<Type> &arg2, const vector<Type> &arg3, int arg4)	\
{												\
	int n2 = arg2.size();									\
	int n3 = arg3.size();									\
	int n = fmax(n2,n3);									\
	vector<Type> res(n);									\
	for(int i=0;i<n;i++) res[i] = FUN(arg1,arg2[i],arg3[i],arg4);				\
	return res;										\
}

/**	\brief Vectorization macro 4.TtTi
	
	For four-arguments functions (Type, Type, Type, int). Vectorize first and third arguments.
	*/
#define VECTORIZE4_TtTi(FUN)									\
template <class Type>										\
vector<Type> FUN(const vector<Type> &arg1, Type arg2, const vector<Type> &arg3, int arg4)	\
{												\
	int n1 = arg1.size();									\
	int n3 = arg3.size();									\
	int n = fmax(n1,n3);									\
	vector<Type> res(n);									\
	for(int i=0;i<n;i++) res[i] = FUN(arg1[i],arg2,arg3[i],arg4);				\
	return res;										\
}

/**	\brief Vectorization macro 4.TTti
	
	For four-arguments functions (Type, Type, Type, int). Vectorize first and second arguments.
	*/
#define VECTORIZE4_TTti(FUN)									\
template <class Type>										\
vector<Type> FUN(const vector<Type> &arg1, const vector<Type> &arg2, Type arg3, int arg4)	\
{												\
	int n1 = arg1.size();									\
	int n2 = arg2.size();									\
	int n = fmax(n1,n2);									\
	vector<Type> res(n);									\
	for(int i=0;i<n;i++) res[i] = FUN(arg1[i],arg2[i],arg3,arg4);				\
	return res;										\
}

/**	\brief Vectorization macro 4.TTTi
	
	For four-arguments functions (Type, Type, Type, int). Vectorize first, second and third arguments.
	*/
#define VECTORIZE4_TTTi(FUN)											\
template <class Type>												\
vector<Type> FUN(const vector<Type> &arg1, const vector<Type> &arg2, const vector<Type> &arg3, int arg4)	\
{														\
	int n1 = arg1.size();											\
	int n2 = arg2.size();											\
	int n3 = arg3.size();											\
	int n = fmax(fmax(n1,n2),n3);										\
	vector<Type> res(n);for(int i=0;i<n;i++) res[i] = FUN(arg1[i],arg2[i],arg3[i],arg4);			\
	return res;												\
}

using CppAD::abs;
VECTORIZE1(abs);
VECTORIZE1(acos);
VECTORIZE1(asin);
VECTORIZE1(atan);
VECTORIZE1(cos);
VECTORIZE1(erf);
VECTORIZE1(exp);
VECTORIZE1(lgamma);
VECTORIZE1(log);
VECTORIZE1(log10);
VECTORIZE1(sin);
VECTORIZE1(sqrt);   //mangler atan2, pow
VECTORIZE2(pow);

VECTORIZE31(dnorm);
VECTORIZE31(dnbinom);
VECTORIZE31(dgamma);
VECTORIZE31(dlgamma);
VECTORIZE21(dpois);


/* max and min of vector */
double max(const vector<double> &x)
{
  double res=x[0];
  for(int i=0;i<x.size();i++){
    //res=CondExpGt(res,x[i],res,x[i]);
    if(res<x[i])res=x[i];
  }
  return res;
}
template <class Type>
Type max(const vector<Type> &x)
{
  Type res=x[0];
  for(int i=0;i<x.size();i++){
    res=CondExpGt(res,x[i],res,x[i]);
  }
  return res;
}
double min(const vector<double> &x)
{
  double res=x[0];
  for(int i=0;i<x.size();i++){
    //res=CondExpGt(res,x[i],res,x[i]);
    if(res>x[i])res=x[i];
  }
  return res;
}
template <class Type>
Type min(const vector<Type> &x)
{
  Type res=x[0];
  for(int i=0;i<x.size();i++){
    res=CondExpLt(res,x[i],res,x[i]);
  }
  return res;
}


/* Normalize rows and colums in a matrix */
// template<class Type>
// void normalizeRows(matrix<Type> &x){
//   for(int i=0;i<x.rows();i++){
//     matrix_row<matrix<Type> > xr(x,i);
//     xr/=sum(xr);
//   }
// }
// template<class Type>
// void normalizeCols(matrix<Type> &x){
//   for(int i=0;i<x.cols();i++){
//     matrix_column<matrix<Type> > xc(x,i);
//     xc/=sum(xc);
//   }
// }


/* Tape the choice of a matrix-row  */
template<class Type>
vector<Type> select_row(const matrix<Type> &x, Type index)
{
   int n1=x.rows();
   int n2=x.cols();
   vector<Type> res(n2);
   //res.clear();
   for(int j=0;j<n2;j++)res(j)=Type(0);
   for(int i=0;i<n1;i++)
     for(int j=0;j<n2;j++)
        res(j)+=CppAD::CondExpEq(index,Type(i),x(i,j),Type(0));
   return res;
}

/* Tape the choice of a vector-element  */
template<class Type>
Type select_elt(const vector<Type> &x, Type i)
{
   int n=x.size();
   Type res(0.0);
   for(int j=0;j<n;j++)
     {
        res+=CppAD::CondExpEq(i,Type(j),x[j],Type(0));
     }
   return res;
}

/* Generate piecewice constant function-object based on n discontinuity-points
   and n+1 function values */
template<class Type>
class piecewice{
private:
  Type y0;
  vector<Type> x,dy;
  int n;
public:
  bool leftcontinuous;
  piecewice(const vector<Type> &x_, const vector<Type> &y,
	    bool leftcontinuous_=true):
    x(x_),leftcontinuous(leftcontinuous_){
    n=x.size();
    if(n!=y.size()-1)error("piecewice: y must be one longer than x");
    dy.resize(n);
    y0=y(0);
    for(int i=0;i<n;i++)dy(i)=y(i+1)-y(i);
  }
  inline Type operator()(const Type &t){
    Type res=y0;
    if(leftcontinuous){
      for(int i=0;i<n;i++)res+=CppAD::CondExpGt(t,x(i),dy(i),Type(0));
    } else {
      for(int i=0;i<n;i++)res+=CppAD::CondExpGe(t,x(i),dy(i),Type(0));
    }
    return res;
  }
}; 
