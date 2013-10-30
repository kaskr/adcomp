/*
  //////////////
  // Remarks: //
  //////////////

  1. Only implemented for 0-order Forward and 1st order Reverse on nested AD types
     up to AD<AD<double> >. 
  2. Hack: All tapes are recorded with parameter vector of ones! 
     Reason: often zero is invalid (gives nan) e.g. "dnorm" with "sd=0".
     But, "one" could also fail!
  3. The implementation of "rev_jac_sparse" is a hack to avoid the atomic functions
     are "optimized out" when calling ADFun.optimize().

 */


# include "cppad/cppad.hpp"
template<class T>
void printVec(T x){
  CppAD::vector<double> y(x.size());
  for(int i=0;i<x.size();i++)y[i]=asDouble(x[i]);
  std::cout << y ; std::cout.flush();
}

/* 
   Auto generated template						
   "double"-version triggers the tape generator!!!			
 */									


/*
  ======== New idea
  To simplify the notation use abbreviations like
  AD3 = AD<AD<AD<double>>> 
  and assume the name of our function is function_name="dnorm".
  
  At highest nested level (AD3) we evaluate dnorm(AD3), and the usual
  atomic template would ensure that our tape contains litterally
      "dnorm,dnorm,dnorm,... etc"
  The new idea is to modify ATOMIC_TEMPLATE so that pf->Forward(0,) and
  pf->Reverse(1,) are replaced by literally dnormFor (or maybe just dnorm)
  and dnormRev, where dnormFor and dnormRev are _atomic_. That requires that
  we already defined FRONTENDS for dnormFor(AD2) and dnormRev(AD2) at 
  global scope. These frontends should act by simply writing literally 
  "dnorm" and "dnormRev" to the tape. And this could be achieved by having
  one more CPPAD_USER_ATOMIC declaration inside ATOMIC_TEMPLATE in the AD2
  case: Not just a "fun" but also a "funRev".
  
  Only, when reaching the lowest nested level (Base=double), we can no
  longer just refer to the global scope, but we have to actually call
  pf->Forward, pf->Reverse(1), pf->Reverse(2).
  
  To summarize we need

  FRONTENDS:
  AD1   dnorm,  dnormRev,  dnormRevRev
  AD2   dnorm,  dnormRev
  AD3   dnorm

  I think that each instance (e.g. AD2-dnormRev) need to be atomic
  *and* be taped. For instance, to obtain AD1-dnormRevRev the tape AD2-dnormRev
  would be needed.  
  

  Or maybe keep the current ATOMIC_TEMPLATE as it is. 

  What we have now:
  ================
  ATOMIC_TEMPLATE that for any vector->vector function does:
    1. Tapes the code of the function
    2. Associates the tape with an atomic function placed internally in CppAD.

  This should be easy:
  ====================
  1. Generate all the following using the old ATOMIC_TEMPLATE + ATOMIC_FRONTEND:
  AD1   dnorm,  dnormRev,  dnormRevRev
  AD2   dnorm,  dnormRev
  AD3   dnorm
  2. Keep only:
  AD1   dnorm,  dnormRev,  dnormRevRev
  3. Generate all the following using a new ATOMIC_TEMPLATE_LAZY where Forward and
  Reverse calls are replaced by the now well defined atomic functions.
  AD1   dnorm,  dnormRev,  dnormRevRev
  AD2   dnorm_lazy,  dnormRev_lazy
  AD3   dnorm_lazy
 
  Now we may ask:
  ================
  Why not merge the ATOMIC_TEMPLATE and ATOMIC_TEMPLATE_LAZY ?
  I dont think we can because of the direction - look into it!


  Idea for modification
  ====================
  1. Replace pfloc->Forward and pfloc->Reverse by their atomic versions
     "dnorm" and "dnormRev"?


 */

#ifdef _OPENMP
#define NTHREADS omp_get_max_threads()
#define THREADNUM omp_get_thread_num()
#else
#define NTHREADS 1
#define THREADNUM 0
#endif

#define ATOMIC_TEMPLATE(function_name,Base,ns_name,code)		\
namespace ns_name{							\
  using CppAD::vector;							\
  using CppAD::AD;							\
  code									\
  CppAD::ADFun<Base >* pf=NULL;						\
  CppAD::vector< CppAD::ADFun<Base >* > vpf; \
  int n,m;								\
  bool forward(								\
	       size_t                   id ,				\
	       size_t                order ,				\
	       size_t                    n ,				\
	       size_t                    m ,				\
	       const vector<bool>&      vx ,				\
	       vector<bool>&           vzy ,				\
	       const vector<Base >&     tx ,				\
	       vector<Base >&          tzy ){				\
    CppAD::ADFun<Base >* pfloc=vpf[THREADNUM];					\
    assert( tx.size() >= (order+1) * n );				\
    assert( tzy.size() >= (order+1) * m );				\
    size_t j,k;								\
    if( vx.size() > 0 ){						\
      assert( vx.size() >= n );						\
      assert( vzy.size() >= m );					\
      for(k=0;k<m;k++)vzy[k]=true;					\
    }									\
    vector<Base > xtmp(n);						\
    vector<Base > ytmp(m);						\
    for(j=0;j<=order;j++){						\
      for(k=0;k<n;k++)xtmp[k]=tx[k*(order+1)+j];			\
      ytmp=pfloc->Forward(j,xtmp);					\
      for(k=0;k<m;k++)tzy[k*(order+1)+j]=ytmp[k];			\
    }									\
    return true;							\
  }									\
  bool reverse(								\
	       size_t                   id ,				\
	       size_t                order ,				\
	       size_t                    n ,				\
	       size_t                    m ,				\
	       const vector<Base >&     tx ,				\
	       const vector<Base >&    tzy ,				\
	       vector<Base >&           px ,				\
	       const vector<Base >&    pzy){				\
    CppAD::ADFun<Base >* pfloc=vpf[THREADNUM];					\
    assert( tx.size() >= (order+1) * n );				\
    assert( tzy.size() >= (order+1) * m );				\
    assert( px.size() >= (order+1) * n );				\
    assert( pzy.size() >= (order+1) * m );				\
    size_t j, k;							\
    vector<Base > xtmp(n);						\
    vector<Base > ytmp(m);						\
    for(j=0;j<=order;j++){						\
      for(k=0;k<n;k++)xtmp[k]=tx[k*(order+1)+j];			\
      ytmp=pfloc->Forward(j,xtmp);					\
    }									\
    vector<Base > pxtmp(n*(order+1));					\
    pxtmp=pfloc->Reverse(order+1,pzy);					\
    for(k=0;k<pxtmp.size();k++)px[k]=pxtmp[k];				\
    return true;							\
  }									\
  bool for_jac_sparse(							\
		      size_t                               id ,		\
		      size_t                                n ,		\
		      size_t                                m ,		\
		      size_t                                q ,		\
		      const vector< std::set<size_t> >&     r ,		\
		      vector< std::set<size_t> >&           s ) {	\
    return false;							\
  }									\
  bool rev_jac_sparse(							\
		      size_t                               id ,		\
		      size_t                                n ,		\
		      size_t                                m ,		\
		      size_t                                q ,		\
		      vector< std::set<size_t> >&           r ,		\
		      const vector< std::set<size_t> >&     s ) {	\
    for(int i=0;i<r.size();i++)r[i].insert(0);   \
    return true;							\
  }									\
  bool rev_hes_sparse(							\
		      size_t                               id ,		\
		      size_t                                n ,		\
		      size_t                                m ,		\
		      size_t                                q ,		\
		      const vector< std::set<size_t> >&     r ,		\
		      const vector<bool>&                   s ,		\
		      vector<bool>&                         t ,		\
		      const vector< std::set<size_t> >&     u ,		\
		      vector< std::set<size_t> >&           v ) {	\
    return false;							\
  }									\
  CPPAD_USER_ATOMIC(							\
  		    fun               ,					\
  		    CppAD::vector ,					\
  		    Base              ,					\
  		    forward           ,					\
  		    reverse           ,					\
  		    for_jac_sparse    ,					\
  		    rev_jac_sparse    ,					\
  		    rev_hes_sparse   )					\
  void generate_tape(int n_, int m_){					\
  std::cout << "Generating tape\n";     					\
    n=n_; m=m_;								\
    vpf.resize(NTHREADS); \
    CppAD::vector< AD<Base > > x(n);				\
    for(int i=0;i<n;i++)x[i]=AD<Base>(1.0); \
    CppAD::Independent(x);						\
    CppAD::vector< AD<Base > > y(m);				\
    y=function_name(x);							\
    pf= new CppAD::ADFun<Base >(x,y);					\
    pf->optimize(); \
    vpf[0]=pf; \
    for(int i=1;i<NTHREADS;i++){ \
      vpf[i]=new CppAD::ADFun<Base >(); \
      vpf[i]->operator=(*pf); \
    } \
    fun(0,x,y); \
  }									\
}

//#define PRINT_INFO(ns_name) std::cout  << #ns_name << " order: "  << order << " n " << n << " m " << m << "\n"; std::cout.flush();
#define PRINT_INFO(ns_name) 
//#define PRINT_INFO(ns_name) std::cout << #ns_name << omp_get_thread_num() << "\n"; std::cout.flush();

//////////////////////////////////////////////////////////////////////////////////////////
// Atomic template that replace all forward and reverse calls by atomic taped functions //
//////////////////////////////////////////////////////////////////////////////////////////
#define ATOMIC_TEMPLATE_LAZY(Base,ns_name,nsFor,nsForRev)		\
namespace ns_name{							\
  using CppAD::vector;							\
  using CppAD::AD;							\
  int n,m;								\
  bool forward(								\
	       size_t                   id ,				\
	       size_t                order ,				\
	       size_t                    n ,				\
	       size_t                    m ,				\
	       const vector<bool>&      vx ,				\
	       vector<bool>&           vzy ,				\
	       const vector<Base >&     tx ,				\
	       vector<Base >&          tzy ){				\
    PRINT_INFO(ns_name)							\
    assert( tx.size() >= (order+1) * n );				\
    assert( tzy.size() >= (order+1) * m );				\
    size_t j,k;								\
    if( vx.size() > 0 ){						\
      assert( vx.size() >= n );						\
      assert( vzy.size() >= m );					\
      for(k=0;k<m;k++)vzy[k]=true;					\
    }									\
    vector<Base > xtmp(n);						\
    vector<Base > ytmp(m);						\
    for(j=0;j<=order;j++){						\
      for(k=0;k<n;k++)xtmp[k]=tx[k*(order+1)+j];			\
      nsFor::fun(0,xtmp,ytmp);						\
      for(k=0;k<m;k++)tzy[k*(order+1)+j]=ytmp[k];			\
    }									\
    return true;							\
  }									\
  bool reverse(								\
	       size_t                   id ,				\
	       size_t                order ,				\
	       size_t                    n ,				\
	       size_t                    m ,				\
	       const vector<Base >&     tx ,				\
	       const vector<Base >&    tzy ,				\
	       vector<Base >&           px ,				\
	       const vector<Base >&    pzy){				\
    PRINT_INFO(ns_name)							\
    assert( tx.size() >= (order+1) * n );				\
    assert( tzy.size() >= (order+1) * m );				\
    assert( px.size() >= (order+1) * n );				\
    assert( pzy.size() >= (order+1) * m );				\
    size_t j, k;							\
    vector<Base > xtmp(n+m);						\
    for(k=0;k<n;k++)xtmp[k]=tx[k];      \
    for(k=0;k<m;k++)xtmp[k+n]=pzy[k];   \
    vector<Base > pxtmp(n*(order+1));					\
    nsForRev::fun(0,xtmp,pxtmp);					\
    for(k=0;k<pxtmp.size();k++)px[k]=pxtmp[k]; \
    return true;							\
  }									\
  bool for_jac_sparse(							\
		      size_t                               id ,		\
		      size_t                                n ,		\
		      size_t                                m ,		\
		      size_t                                q ,		\
		      const vector< std::set<size_t> >&     r ,		\
		      vector< std::set<size_t> >&           s ) {	\
    return false;							\
  }									\
  bool rev_jac_sparse(							\
		      size_t                               id ,		\
		      size_t                                n ,		\
		      size_t                                m ,		\
		      size_t                                q ,		\
		      vector< std::set<size_t> >&           r ,		\
		      const vector< std::set<size_t> >&     s ) {	\
    for(int i=0;i<r.size();i++)r[i].insert(0);   \
    return true;							\
  }									\
  bool rev_hes_sparse(							\
		      size_t                               id ,		\
		      size_t                                n ,		\
		      size_t                                m ,		\
		      size_t                                q ,		\
		      const vector< std::set<size_t> >&     r ,		\
		      const vector<bool>&                   s ,		\
		      vector<bool>&                         t ,		\
		      const vector< std::set<size_t> >&     u ,		\
		      vector< std::set<size_t> >&           v ) {	\
    return false;							\
  }									\
  CPPAD_USER_ATOMIC(							\
  		    fun               ,					\
  		    CppAD::vector ,					\
  		    Base              ,					\
  		    forward           ,					\
  		    reverse           ,					\
  		    for_jac_sparse    ,					\
  		    rev_jac_sparse    ,					\
  		    rev_hes_sparse   )					\
  void generate_tape(int n_, int m_){					\
  std::cout << "trigger lazy atomic\n";     					\
    n=n_; m=m_;								\
    CppAD::vector< AD<Base > > x(n);				\
    for(int i=0;i<n;i++)x[i]=AD<Base>(1.0); \
    CppAD::vector< AD<Base > > y(m);				\
    fun(0,x,y); \
  }									\
}


/* This is how the atomic function appears at user level,
   for instance if function_name="dnorm" then we get versions
   1. vector<double> dnorm(vector<double>)
   2. vector<AD<double>> dnorm(vector<AD<double>>)
   etc.
 */

// #define ATOMIC_FRONTEND(function_name,Base,ns_name)	\
//   CppAD::vector<CppAD::AD<Base > > function_name	\
//   (CppAD::vector<CppAD::AD<Base > > x){			\
//   assert( NULL!=ns_name::pf );				\
//   assert( x.size()==ns_name::n );			\
//   CppAD::vector<CppAD::AD<Base > > y(ns_name::m);	\
//   ns_name::fun(0,x,y);					\
//   return y;						\
// }									


/* Frontend + wrapper tmbutils::vector->CppAD::vector*/
#define ATOMIC_FRONTEND(function_name,Base,ns_name)	\
  CppAD::vector<CppAD::AD<Base > > function_name	\
  (CppAD::vector<CppAD::AD<Base > > x){			\
  CppAD::vector<CppAD::AD<Base > > y(ns_name::m);	\
  ns_name::fun(0,x,y);					\
  return y;						\
}							\
tmbutils::vector<CppAD::AD<Base > > function_name		\
  (tmbutils::vector<CppAD::AD<Base > > x){			\
    return tmbutils::vector<CppAD::AD<Base > >(function_name(	\
    CppAD::vector<CppAD::AD<Base > >(x)));		\
}									



/* Trigger tape-generation + wrapper tmbutils::vector->CppAD::vector*/
#define AD1 CppAD::AD<double>
#define AD2 CppAD::AD<CppAD::AD<double> >
#define ATOMIC_TRIGGER(function_name,ns_name,ns2,ns3)			\
CppAD::vector<double> function_name(CppAD::vector<double> x){		\
  CppAD::vector<double> y=ns_name::function_name(x);			\
  if(ns_name::pf==NULL){						\
    int n=x.size();							\
    int m=y.size();							\
    ns_name::generate_tape(n,m);					\
    ns2::generate_tape(n,m);						\
    ns3::generate_tape(n,m);						\
    atomic4::generate_tape(m+n,n);					\
    atomic5::generate_tape(m+n,n);					\
    atomic6::generate_tape(m+n+n,m+n);					\
    lazy2::generate_tape(n,m);						\
    lazy3::generate_tape(n,m);						\
    lazy5::generate_tape(m+n,n);					\
  }									\
  return y;								\
}									\
tmbutils::vector<double > function_name					\
(tmbutils::vector<double > x){						\
  return tmbutils::vector<double >(function_name(CppAD::vector<double >(x))); \
}									




//template<class T>T myf(T x){return atomic2::pf->Reverse(1,x);}
/*
  Full taped version:
  ===================
  AD1   dnorm,   dnormRev,  dnormRevRev
  AD2   dnorm,   dnormRev
  AD3   dnorm
  ...located in these namespaces:
        atomic,  atomic4,   atomic6
	atomic2, atomic5
	atomic3,

  Lazy versions:
  ==============
  ...located in these namespaces:
        lazy,    lazy4,     lazy6   (<-- Equal to full taped versions)
	lazy2,   lazy5              (<-- E.g. lazy2 calls lazy and lazy4)
	lazy3,                      (<-- E.g. lazy3 calls lazy2 and lazy5)

*/
#define FORREW(ns)				\
template<class T>T myf(T x){			\
  int n=ns::n; int m=ns::m;			\
  T tmp1(n);					\
  T tmp2(m);					\
  for(size_t i=0;i<n;i++)tmp1[i]=x[i];		\
  for(size_t i=n;i<n+m;i++)tmp2[i-n]=x[i];	\
  ns::vpf[0]->Forward(0,tmp1);			\
  return ns::vpf[0]->Reverse(1,tmp2);		\
}

#define ATOMIC_FUNCTION(function_name,code)				\
ATOMIC_TEMPLATE(function_name,double,atomic,code)			\
ATOMIC_TEMPLATE(function_name,CppAD::AD<double>,atomic2,code)		\
ATOMIC_TEMPLATE(function_name,CppAD::AD<CppAD::AD<double> >,atomic3,code) \
ATOMIC_TEMPLATE(myf,double,atomic4,FORREW(atomic2))			\
ATOMIC_TEMPLATE(myf,AD1,atomic5,FORREW(atomic3))			\
ATOMIC_TEMPLATE(myf,double,atomic6,FORREW(atomic5))			\
ATOMIC_TEMPLATE_LAZY(AD1,lazy2,atomic,atomic4)				\
ATOMIC_TEMPLATE_LAZY(AD1,lazy5,atomic4,atomic6)				\
ATOMIC_TEMPLATE_LAZY(AD2,lazy3,lazy2,lazy5)				\
ATOMIC_FRONTEND(function_name,double ,atomic)				\
ATOMIC_FRONTEND(function_name,CppAD::AD<double>,lazy2)			\
ATOMIC_FRONTEND(function_name,CppAD::AD<CppAD::AD<double> >,lazy3)	\
ATOMIC_TRIGGER(function_name,atomic,atomic2,atomic3)



/* And finally, all the user needs:

   1. Define a template function mapping vector to vector, e.g
   template<class T>
   vector<T> foo(vector<T>);

   2. REGISTER_ATOMIC(foo);

*/
#define REGISTER_ATOMIC(function_name)				\
ATOMIC_FUNCTION(function_name,					\
template<class Type>						\
CppAD::vector<Type> function_name(CppAD::vector<Type> x){	\
  return ::function_name(tmbutils::vector<Type>(x));			\
})

