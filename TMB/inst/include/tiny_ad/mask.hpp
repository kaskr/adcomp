template<int base, int n>
struct pow_t {
  static const long int value = base * pow_t<base, n-1>::value;
};
template<int base>
struct pow_t<base, 0> {
  static const long int value = 1;
};
template<long int mask>
struct mask_t {
  template<int length, int i=0>
  struct set_length {
    static const int base = 8;
    static const long int power = pow_t<base, i>::value;
    typedef set_length<length, i+1> next_index_t;
    static const int flag = ( (mask / power) % base ) != 0;
    static const int count = flag + next_index_t::count;
    static const int Id = count - 1;
    static const int Index = length - i - 1;
    next_index_t next_index;
    template<class S, class T>
    void copy(S &dest, const T &orig) {
      dest[Index] = (flag ? orig[Id] : 0);
      next_index.copy(dest, orig);
    }
    template<class S, class T>
    void activate_derivs(S &var, T &value) {
      var[Index] = value[Index];
      if (flag) var[Index].setid(Id);
      next_index.activate_derivs(var, value);
    }
  };
  template<int length>
  struct set_length<length, length> {
    static const int count = 0;
    void trace() { }
    template<class S, class T>
    void copy(S &dest, const T &orig) { }
    template<class S, class T>
    void activate_derivs(S &var, T &value) { }
  };
};

template<int nvar>
struct ADTypes {
  typedef tiny_ad::variable<1, nvar> order1;
  typedef tiny_ad::variable<2, nvar> order2;
  typedef tiny_ad::variable<3, nvar> order3;
};

// 'TMB_BIND_ATOMIC' depends on these:
#define NCHAR(x) sizeof(#x)-1
#define OCTAL(x) 0 ## x

/** \brief Bind an atomic function to a forward differentiated
    template function using tiny_ad. Generates up to 3rd order
    derivatives.

    \param NAME Name of atomic function.

    \param MASK Binary mask denoting active variables. e.g. 0101 says
    there are a total of four input parameters of which two are
    active.

    \param CALL A call to a differentiable template function. The
    input parameters must be referred to as 'x'.

    \note: The argument list is expanded with a number denoting the
    order. So, in a case with four input parameters the generated
    atomic function actually has five parameters. The last parameter
    must be set to zero when calling the atomic function.
*/
#define TMB_BIND_ATOMIC(NAME,MASK,CALL)					\
TMB_ATOMIC_VECTOR_FUNCTION(						\
  NAME									\
  ,									\
  (size_t)     								\
  pow((double)								\
  atomic::mask_t<OCTAL(MASK)>::set_length< NCHAR(MASK) >::count,	\
      CppAD::Integer(tx[NCHAR(MASK)]))					\
  ,									\
  int order = CppAD::Integer(tx[NCHAR(MASK)]);				\
  typedef								\
  atomic::mask_t<OCTAL(MASK)>::set_length<NCHAR(MASK)> mask_type;	\
  mask_type mask;							\
  static const int nvar = mask_type::count;				\
  atomic::tiny_vec_ref<double> tyref(&ty[0], ty.size());                \
  if(order==0) {							\
    typedef double Float;						\
    CppAD::vector<Float> x(tx);						\
    ty[0] = CALL;							\
  }									\
  else if (order==1) {							\
    typedef typename atomic::ADTypes<nvar>::order1 Float;		\
    Float x[NCHAR(MASK)];						\
    mask.activate_derivs(x, tx);					\
    tyref = CALL.getDeriv();                                            \
  }									\
  else if (order==2) {							\
    typedef typename atomic::ADTypes<nvar>::order2 Float;		\
    Float x[NCHAR(MASK)];						\
    mask.activate_derivs(x, tx);					\
    tyref = CALL.getDeriv();                                            \
  }									\
  else if (order==3) {							\
    typedef typename atomic::ADTypes<nvar>::order3 Float;		\
    Float x[NCHAR(MASK)];						\
    mask.activate_derivs(x, tx);					\
    tyref = CALL.getDeriv();                                            \
  }									\
  else									\
    Rf_error("Order not implemented");					\
  ,									\
  typedef								\
  atomic::mask_t<OCTAL(MASK)>::set_length<NCHAR(MASK)> mask_type;	\
  mask_type mask;							\
  static const int nvar = mask_type::count;				\
  CppAD::vector<Type> tx_(tx);						\
  tx_[NCHAR(MASK)] = tx_[NCHAR(MASK)] + Type(1.0);			\
  vector<Type> tmp = NAME(tx_);						\
  matrix<Type> m = tmp.matrix();					\
  m.resize(nvar, m.size() / nvar);					\
  vector<Type> w = py;							\
  vector<Type> px_ = m * w.matrix();					\
  mask.copy(px, px_);							\
  px[NCHAR(MASK)] = 0;							\
  )

// ======================================================================


#ifdef TMBAD_FRAMEWORK

#undef TMB_BIND_ATOMIC
#ifndef TMB_MAX_ORDER
#define TMB_MAX_ORDER 3
#endif

#define TMB_BIND_ATOMIC(NAME,MASK,CALL)                                 \
template<int order, int ninput, int noutput, long int mask>             \
struct NAME ## Eval {                                                   \
  typedef typename                                                      \
  atomic::mask_t<mask>::template set_length<ninput> mask_type;          \
  mask_type mask_;                                                      \
  static const int nvar = mask_type::count;                             \
  template <class S, class T>                                           \
  void operator()(S* tx, T* ty) {                                       \
    typedef atomic::tiny_ad::variable<order, nvar> Float;               \
    atomic::tiny_vec_ref<double> tyref(&(ty[0]), noutput);              \
    Float x[ninput];                                                    \
    mask_.activate_derivs(x, tx);                                       \
    tyref = (CALL).getDeriv();                                          \
  }                                                                     \
};                                                                      \
template<int ninput, int noutput, long int mask>                        \
struct NAME ## Eval<0, ninput, noutput, mask> {                         \
  template <class S, class T>                                           \
  void operator()(S* tx, T* ty) {                                       \
    S* x = tx;                                                          \
    ty[0] = (CALL);                                                     \
  }                                                                     \
};                                                                      \
template<int order, int ninput, int noutput, long int mask>             \
struct NAME ## Op : TMBad::global::Operator<ninput, noutput> {          \
  static const bool add_forward_replay_copy = true;                     \
  typedef typename                                                      \
  atomic::mask_t<mask>::template set_length<ninput> mask_type;          \
  mask_type mask_;                                                      \
  static const int nvar = mask_type::count;                             \
  template <class S, class T>                                           \
  void eval(S* tx, T* ty) const {                                       \
    NAME ## Eval<order, ninput, noutput, mask>()(tx, ty);               \
  }                                                                     \
  std::vector<TMBad::ad_plain>                                          \
  add_to_tape(const std::vector<TMBad::ad_plain> &x) {                  \
    TMBad::OperatorPure* pOp = TMBad::get_glob()->getOperator<NAME ## Op>(); \
    return                                                              \
      TMBad::get_glob()->add_to_stack<NAME ## Op>(pOp, x);              \
  }                                                                     \
  std::vector<TMBad::ad_plain>                                          \
  operator()(const std::vector<TMBad::ad_plain> &x) {                   \
    return add_to_tape(x);                                              \
  }                                                                     \
  Eigen::Matrix<TMBad::ad_aug, nvar, noutput / nvar>                    \
  operator()(const Eigen::Array<TMBad::ad_aug, ninput, 1> &x) {         \
    std::vector<TMBad::ad_plain> x_(&(x(0)), &(x(0)) + x.size());       \
    Eigen::Matrix<TMBad::ad_aug, nvar, noutput / nvar> ans;             \
    std::vector<TMBad::ad_plain> y = add_to_tape(x_);                   \
    for (size_t i=0; i<y.size(); i++) ans(i) = y[i];                    \
    return ans;                                                         \
  }                                                                     \
  Eigen::Matrix<double, nvar, noutput / nvar>                           \
  operator()(const Eigen::Array<double, ninput, 1> &x) {                \
    Eigen::Matrix<double, nvar, noutput / nvar> ans;                    \
    eval(&(x(0)), &(ans(0)));                                           \
    return ans;                                                         \
  }                                                                     \
  template<class Type>                                                  \
  void forward(TMBad::ForwardArgs<Type> &args) {                        \
    Rf_error("Un-implemented method request");                          \
  }                                                                     \
  void forward(TMBad::ForwardArgs<double> &args) {                      \
    double x[ninput];                                                   \
    for (size_t i=0; i<ninput; i++) x[i] = args.x(i);                   \
    eval(x, &(args.y(0)));                                              \
  }                                                                     \
  template<class Type>                                                  \
  void reverse(TMBad::ReverseArgs<Type> &args) {                        \
    Eigen::Array<Type, ninput, 1> tx;                                   \
    for (size_t i=0; i<ninput; i++) tx(i) = args.x(i);                  \
    Eigen::Matrix<Type, noutput, 1> w;                                  \
    for (size_t i=0; i<noutput; i++) w(i) = args.dy(i);                 \
    NAME ## Op<order+1, ninput, noutput * nvar, mask> foo;              \
    Eigen::Matrix<Type, nvar, noutput> ty;                              \
    ty = foo(tx);                                                       \
    Eigen::Matrix<Type, nvar, 1> tyw = ty * w;                          \
    Type tmp[ninput];                                                   \
    mask_.copy(tmp, &(tyw[0]));                                         \
    for (size_t i=0; i<ninput; i++) args.dx(i) += tmp[i];               \
  }                                                                     \
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) {               \
    Rf_error("Un-implemented method request");                          \
  }                                                                     \
  const char* op_name() { return #NAME ; }                              \
};                                                                      \
template<int ninput, int noutput, long int mask>                        \
struct NAME ## Op<TMB_MAX_ORDER+1, ninput, noutput, mask> {             \
  typedef typename                                                      \
  atomic::mask_t<mask>::template set_length<ninput> mask_type;          \
  mask_type mask_;                                                      \
  static const int nvar = mask_type::count;                             \
  Eigen::Matrix<TMBad::ad_aug, nvar, noutput / nvar>                    \
  operator()(const Eigen::Array<TMBad::ad_aug, ninput, 1> &x) {         \
    Eigen::Matrix<TMBad::ad_aug, nvar, noutput / nvar> ans;             \
    Rf_error("Order not implemented. Please increase TMB_MAX_ORDER");   \
    return ans;                                                         \
  }                                                                     \
  Eigen::Matrix<double, nvar, noutput / nvar>                           \
  operator()(const Eigen::Array<double, ninput, 1> &x) {                \
    Eigen::Matrix<double, nvar, noutput / nvar> ans;                    \
    Rf_error("Order not implemented. Please increase TMB_MAX_ORDER");   \
    return ans;                                                         \
  }                                                                     \
};                                                                      \
template<class dummy=void>                                              \
CppAD::vector<double>                                                   \
NAME (const CppAD::vector<double> &x) CSKIP_ATOMIC({                    \
  int n = x.size() - 1;                                                 \
  int order = CppAD::Integer(x[n]);                                     \
  typedef NAME ## Op<0, NCHAR(MASK),    1, OCTAL(MASK)> Foo0;           \
  static const int nvar = Foo0::nvar;                                   \
  typedef NAME ## Op<1, NCHAR(MASK), nvar, OCTAL(MASK)> Foo1;           \
  if (order==0) {                                                       \
    CppAD::vector<double> y(1);                                         \
    y[0] = CALL;                                                        \
    return y;                                                           \
  }                                                                     \
  else if (order==1) {                                                  \
    Foo1 foo1;                                                          \
    CppAD::vector<double> y(nvar);                                      \
    foo1.eval(&x[0], &y[0]);                                            \
    return y;                                                           \
  }                                                                     \
  else {                                                                \
    Rf_error("This interface is limited to 0th and 1st deriv order");   \
  }                                                                     \
})                                                                      \
template<class dummy=void>                                              \
CppAD::vector<TMBad::ad_aug>                                            \
NAME (const CppAD::vector<TMBad::ad_aug> &x) CSKIP_ATOMIC({             \
  bool all_constant = true;                                             \
  for (size_t i = 0; i<x.size(); i++)                                   \
    all_constant &= x[i].constant();                                    \
  if (all_constant) {                                                   \
    CppAD::vector<double> xd(x.size());                                 \
    for (size_t i=0; i<xd.size(); i++) xd[i] = x[i].Value();            \
    CppAD::vector<double> yd = NAME(xd);                                \
    CppAD::vector<TMBad::ad_aug> y(yd.size());                          \
    for (size_t i=0; i<yd.size(); i++) y[i] = yd[i];                    \
    return y;                                                           \
  }                                                                     \
  int n = x.size() - 1;                                                 \
  int order = CppAD::Integer(x[n]);                                     \
  std::vector<TMBad::ad_plain> x_(&(x[0]), &(x[0]) + n);                \
  std::vector<TMBad::ad_plain> y_;                                      \
  typedef NAME ## Op<0, NCHAR(MASK),    1, OCTAL(MASK)> Foo0;           \
  static const int nvar = Foo0::nvar;                                   \
  typedef NAME ## Op<1, NCHAR(MASK), nvar, OCTAL(MASK)> Foo1;           \
  if (order==0) {                                                       \
    Foo0 foo0;                                                          \
    y_ = foo0(x_);                                                      \
  }                                                                     \
  else if (order==1) {                                                  \
    Foo1 foo1;                                                          \
    y_ = foo1(x_);                                                      \
  }                                                                     \
  else {                                                                \
    Rf_error("This interface is limited to 0th and 1st deriv order");   \
  }                                                                     \
  CppAD::vector<TMBad::ad_aug> y(y_.size());                            \
  for (size_t i=0; i<y.size(); i++) y[i] = y_[i];                       \
  return y;                                                             \
})                                                                      \
IF_TMB_PRECOMPILE_ATOMICS(                                              \
template                                                                \
CppAD::vector<TMBad::ad_aug>                                            \
NAME<> (const CppAD::vector<TMBad::ad_aug> &x);                         \
template                                                                \
CppAD::vector<double>                                                   \
NAME<> (const CppAD::vector<double> &x);                                \
)

#endif // TMBAD_FRAMEWORK
