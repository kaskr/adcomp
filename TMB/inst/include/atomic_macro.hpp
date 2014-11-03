/* Conditionally skip compilation */
#ifdef WITH_LIBTMB
#define CSKIP(x) ;
#endif
#ifndef WITH_LIBTMB
#define CSKIP(x) x
#endif

#define TMB_ATOMIC_VECTOR_FUNCTION(ATOMIC_NAME,OUTPUT_DIM,ATOMIC_DOUBLE,ATOMIC_REVERSE) \
CppAD::vector<double> ATOMIC_NAME(CppAD::vector<double> tx)CSKIP({	\
  CppAD::vector<double> ty(OUTPUT_DIM);					\
  ATOMIC_DOUBLE;							\
  return ty;								\
})									\
template <class Type>							\
CppAD::vector<AD<Type > > ATOMIC_NAME(CppAD::vector<AD<Type> > x);	\
template <class Type>							\
class atomic##ATOMIC_NAME : public CppAD::atomic_base<Type> {		\
public:									\
  atomic##ATOMIC_NAME(const char* name) : CppAD::atomic_base<Type>(name){ \
    std::cout << "Constructing atomic " << #ATOMIC_NAME << "\n" ;	\
    this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);		\
  }									\
private:								\
  virtual bool forward(size_t p,					\
		       size_t q,					\
		       const CppAD::vector<bool>& vx,			\
		       CppAD::vector<bool>& vy,				\
		       const CppAD::vector<Type>& tx,			\
		       CppAD::vector<Type>& ty				\
		       )						\
  {									\
    if(q>0)error("Atomic '" #ATOMIC_NAME "' order not implemented.\n");	\
    if( vx.size() > 0 ){						\
      bool anyvx = false;						\
      for(int i=0;i<vx.size();i++)anyvx |= vx[i];			\
      for(int i=0;i<vy.size();i++)vy[i] = anyvx;			\
    }									\
    ty = ATOMIC_NAME(tx);	       					\
    return true;							\
  }									\
  virtual bool reverse(size_t q,					\
		       const CppAD::vector<Type>& tx,			\
		       const CppAD::vector<Type>& ty,			\
		       CppAD::vector<Type>& px,				\
		       const CppAD::vector<Type>& py			\
		       )						\
  {									\
    if(q>0)error("Atomic '" #ATOMIC_NAME "' order not implemented.\n");	\
    ATOMIC_REVERSE;							\
    return true;							\
  }									\
  virtual bool rev_sparse_jac(size_t q,					\
			      const CppAD::vector<bool>& rt,		\
			      CppAD::vector<bool>& st)			\
  {									\
    bool anyrt = false;							\
    for(int i=0;i<rt.size();i++)anyrt |= rt[i];				\
    for(size_t i=0;i<st.size();i++)st[i]=anyrt;				\
    return true;							\
  }									\
  virtual bool rev_sparse_jac(size_t q,					\
			      const CppAD::vector< std::set<size_t> >& rt, \
			      CppAD::vector< std::set<size_t> >& st)	\
  {									\
    error("Should not be called");					\
  }									\
};									\
template<class Type> 							\
CppAD::vector<AD<Type > > ATOMIC_NAME(CppAD::vector<AD<Type > > tx){	\
  static atomic##ATOMIC_NAME<Type > afun##ATOMIC_NAME("atomic_" #ATOMIC_NAME); \
  CppAD::vector<AD<Type > > ty(OUTPUT_DIM);				\
  afun##ATOMIC_NAME(tx,ty);						\
  return ty;								\
}


