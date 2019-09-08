
#ifdef TMB_ATOMIC_VECTOR_FUNCTION
#undef TMB_ATOMIC_VECTOR_FUNCTION
#endif

/** \brief Construct atomic vector function based on known derivatives */
#define TMB_ATOMIC_VECTOR_FUNCTION(                                     \
  ATOMIC_NAME, OUTPUT_DIM,                                              \
  ATOMIC_DOUBLE,                                                        \
  ATOMIC_REVERSE                                                        \
)                                                                       \
template<class dummy=void>                                              \
CppAD::vector<TMBad::ad_aug> ATOMIC_NAME (const CppAD::vector<TMBad::ad_aug> &x); \
template<class dummy=void>                                              \
CppAD::vector<double> ATOMIC_NAME (const CppAD::vector<double> &tx) {   \
  CppAD::vector<double> ty(OUTPUT_DIM);                                 \
  ATOMIC_DOUBLE;                                                        \
  return ty;                                                            \
}                                                                       \
template<class dummy=void>                                              \
struct ATOMIC_NAME ## Op : TMBad::global::DynamicInputOutputOperator {  \
  typedef TMBad::global::DynamicInputOutputOperator Base;               \
  ATOMIC_NAME ## Op (TMBad::Index n, TMBad::Index m) : Base(n, m) {}    \
  const char* op_name() { return #ATOMIC_NAME; }                        \
  void forward(TMBad::ForwardArgs<TMBad::Scalar> _args_) {              \
    CppAD::vector<TMBad::Scalar> tx(this->input_size());                \
    CppAD::vector<TMBad::Scalar> ty(this->output_size());               \
    for (size_t i=0; i<tx.size(); i++) tx[i] = _args_.x(i);             \
    ATOMIC_DOUBLE;                                                      \
    for (size_t i=0; i<ty.size(); i++) _args_.y(i) = ty[i];             \
  }                                                                     \
  void forward(TMBad::ForwardArgs<TMBad::Replay> _args_) {              \
    CppAD::vector<TMBad::Replay> tx(this->input_size());                \
    for (size_t i=0; i<tx.size(); i++) tx[i] = _args_.x(i);             \
    CppAD::vector<TMBad::Replay> ty = ATOMIC_NAME(tx);                  \
    for (size_t i=0; i<ty.size(); i++) _args_.y(i) = ty[i];             \
  }                                                                     \
  template<class Type> void reverse(TMBad::ReverseArgs<Type> _args_) {  \
    CppAD::vector<Type> tx(this->input_size());                         \
    CppAD::vector<Type> ty(this->output_size());                        \
    CppAD::vector<Type> px(this->input_size());                         \
    CppAD::vector<Type> py(this->output_size());                        \
    for (size_t i=0; i<tx.size(); i++) tx[i] = _args_.x(i);             \
    for (size_t i=0; i<ty.size(); i++) ty[i] = _args_.y(i);             \
    for (size_t i=0; i<py.size(); i++) py[i] = _args_.dy(i);            \
    ATOMIC_REVERSE;                                                     \
    for (size_t i=0; i<px.size(); i++) _args_.dx(i) += px[i];           \
  }                                                                     \
  void forward(TMBad::ForwardArgs<TMBad::Writer> &args) { ASSERT(false); } \
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) { ASSERT(false); } \
};                                                                      \
template<class dummy=void>                                              \
CppAD::vector<TMBad::ad_aug> ATOMIC_NAME (const CppAD::vector<TMBad::ad_aug> &tx) { \
  TMBad::Index n = tx.size();                                           \
  TMBad::Index m = OUTPUT_DIM;                                          \
  typedef ATOMIC_NAME ## Op <> OP;                                      \
  bool all_constant = true;                                             \
  for (size_t i = 0; i<tx.size(); i++) all_constant &= tx[i].constant(); \
  CppAD::vector<TMBad::ad_aug> ty(m);                                   \
  if (all_constant) {                                                   \
    CppAD::vector<double> xd(tx.size());                                \
    for (size_t i=0; i<xd.size(); i++) xd[i] = tx[i].Value();           \
    CppAD::vector<double> yd = ATOMIC_NAME(xd);                         \
    for (size_t i=0; i<yd.size(); i++) ty[i] = yd[i];                   \
  } else {                                                              \
    TMBad::OperatorPure* pOp = TMBad::get_glob()->getOperator<OP>(n, m); \
    std::vector<TMBad::ad_plain> x(&tx[0], &tx[0] + tx.size());         \
    std::vector<TMBad::ad_plain> y = TMBad::get_glob()->add_to_stack<OP>(pOp, x); \
    for (size_t i=0; i<y.size(); i++) ty[i] = y[i];                     \
  }                                                                     \
  return ty;                                                            \
}                                                                       \
template<class dummy=void>                                              \
void ATOMIC_NAME (const CppAD::vector<TMBad::ad_aug> &tx,               \
                  CppAD::vector<TMBad::ad_aug> &ty) {                   \
  ty = ATOMIC_NAME(tx);                                                 \
}
