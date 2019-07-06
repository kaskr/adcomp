
#ifdef TMB_ATOMIC_VECTOR_FUNCTION
#undef TMB_ATOMIC_VECTOR_FUNCTION
#endif

/** \brief Construct atomic vector function based on known derivatives */
#define TMB_ATOMIC_VECTOR_FUNCTION(                                     \
  ATOMIC_NAME, OUTPUT_DIM,                                              \
  ATOMIC_DOUBLE,                                                        \
  ATOMIC_REVERSE                                                        \
)                                                                       \
CppAD::vector<TMBad::ad_aug> ATOMIC_NAME (const CppAD::vector<TMBad::ad_aug> &x); \
CppAD::vector<double> ATOMIC_NAME (const CppAD::vector<double> &tx) {   \
  CppAD::vector<double> ty(OUTPUT_DIM);                                 \
  ATOMIC_DOUBLE;                                                        \
  return ty;                                                            \
}                                                                       \
struct ATOMIC_NAME ## Op : TMBad::global::DynamicInputOutputOperator {  \
  typedef TMBad::global::DynamicInputOutputOperator Base;               \
  ATOMIC_NAME ## Op (TMBad::Index n, TMBad::Index m) : Base(n, m) {}    \
  const char* op_name() { return "shit"; }                              \
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
    for (size_t i=0; i<py.size(); i++) _args_.dx(i) += py[i];           \
  }                                                                     \
  void forward(TMBad::ForwardArgs<TMBad::Writer> &args) { ASSERT(false); } \
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) { ASSERT(false); } \
};                                                                      \
CppAD::vector<TMBad::ad_aug> ATOMIC_NAME (const CppAD::vector<TMBad::ad_aug> &tx) { \
  TMBad::Index n = tx.size();                                           \
  TMBad::Index m = OUTPUT_DIM;                                          \
  typedef ATOMIC_NAME ## Op OP;                                         \
  TMBad::OperatorPure* pOp = TMBad::get_glob()->getOperator<OP>(n, m);  \
  std::vector<TMBad::ad_plain> x(&tx[0], &tx[0] + tx.size());           \
  std::vector<TMBad::ad_plain> y = TMBad::get_glob()->add_to_stack<OP>(pOp, x); \
  CppAD::vector<TMBad::ad_aug> ty(y.size());                            \
  for (size_t i=0; i<y.size(); i++) ty[i] = y[i];                       \
  return ty;                                                            \
}                                                                       \
void ATOMIC_NAME (const CppAD::vector<TMBad::ad_aug> &tx,               \
                  CppAD::vector<TMBad::ad_aug> &ty) {                   \
  ty = ATOMIC_NAME(tx);                                                 \
}
