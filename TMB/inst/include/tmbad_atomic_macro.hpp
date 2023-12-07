

#define TMB_ATOMIC_VECTOR_FUNCTION_DEFINE(                              \
  ATOMIC_NAME, OUTPUT_DIM,                                              \
  ATOMIC_DOUBLE,                                                        \
  ATOMIC_REVERSE                                                        \
)                                                                       \
template<class dummy>                                                   \
CppAD::vector<TMBad::ad_aug>                                            \
ATOMIC_NAME (const CppAD::vector<TMBad::ad_aug> &x);                    \
template<class dummy>                                                   \
CppAD::vector<double>                                                   \
ATOMIC_NAME (const CppAD::vector<double> &tx) CSKIP_ATOMIC({            \
  CppAD::vector<double> ty(OUTPUT_DIM);                                 \
  ATOMIC_DOUBLE;                                                        \
  return ty;                                                            \
})                                                                      \
template<class dummy=void>                                              \
struct ATOMIC_NAME ## Op : TMBad::global::DynamicInputOutputOperator {  \
  typedef TMBad::global::DynamicInputOutputOperator Base;               \
  ATOMIC_NAME ## Op (TMBad::Index n, TMBad::Index m) : Base(n, m) {}    \
  const char* op_name() { return #ATOMIC_NAME; }                        \
  static const bool add_static_identifier = true;                       \
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
    if (isDouble<Type>::value &&                                        \
          this->output_size() == 1 &&                                   \
            _args_.dy(0) == Type(0)) {                                  \
              return;                                                   \
    }                                                                   \
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
  void forward                                                          \
  (TMBad::ForwardArgs<TMBad::Writer> &args) { TMBAD_ASSERT(false); }    \
  void reverse                                                          \
  (TMBad::ReverseArgs<TMBad::Writer> &args) { TMBAD_ASSERT(false); }    \
};                                                                      \
template<class dummy>                                                   \
CppAD::vector<TMBad::ad_aug>                                            \
ATOMIC_NAME (const CppAD::vector<TMBad::ad_aug> &tx) CSKIP_ATOMIC({     \
  TMBad::Index n = tx.size();                                           \
  TMBad::Index m = OUTPUT_DIM;                                          \
  typedef ATOMIC_NAME ## Op <> OP;                                      \
  bool all_constant = true;                                             \
  for (size_t i = 0; i<tx.size(); i++)                                  \
    all_constant &= tx[i].constant();                                   \
  CppAD::vector<TMBad::ad_aug> ty(m);                                   \
  if (all_constant) {                                                   \
    CppAD::vector<double> xd(tx.size());                                \
    for (size_t i=0; i<xd.size(); i++) xd[i] = tx[i].Value();           \
    CppAD::vector<double> yd = ATOMIC_NAME(xd);                         \
    for (size_t i=0; i<yd.size(); i++) ty[i] = yd[i];                   \
  } else {                                                              \
    TMBad::OperatorPure*                                                \
      pOp = TMBad::get_glob()->getOperator<OP>(n, m);                   \
    std::vector<TMBad::ad_plain>                                        \
      x(&tx[0], &tx[0] + tx.size());                                    \
    std::vector<TMBad::ad_plain>                                        \
      y = TMBad::get_glob()->add_to_stack<OP>(pOp, x);                  \
    for (size_t i=0; i<y.size(); i++) ty[i] = y[i];                     \
  }                                                                     \
  return ty;                                                            \
})                                                                      \
template<class dummy=void>                                              \
void ATOMIC_NAME (const CppAD::vector<TMBad::ad_aug> &tx,               \
                  CppAD::vector<TMBad::ad_aug> &ty) {                   \
  ty = ATOMIC_NAME(tx);                                                 \
}                                                                       \
IF_TMB_PRECOMPILE_ATOMICS(                                              \
template                                                                \
CppAD::vector<double>                                                   \
ATOMIC_NAME<> (const CppAD::vector<double>& tx);                        \
template                                                                \
CppAD::vector<TMBad::ad_aug>                                            \
ATOMIC_NAME<>(const CppAD::vector<TMBad::ad_aug>& tx);                  \
)

#define TMB_ATOMIC_STATIC_FUNCTION(                                     \
  ATOMIC_NAME,                                                          \
  INPUT_SIZE,                                                           \
  ATOMIC_DOUBLE,                                                        \
  ATOMIC_REVERSE                                                        \
)                                                                       \
template<class dummy=void>                                              \
CppAD::vector<TMBad::ad_aug> ATOMIC_NAME                                \
(const CppAD::vector<TMBad::ad_aug> &x);                                \
template<class dummy=void>                                              \
CppAD::vector<double> ATOMIC_NAME                                       \
(const CppAD::vector<double> &tx) CSKIP_ATOMIC({                        \
  CppAD::vector<double> ty(1);                                          \
  ATOMIC_DOUBLE;                                                        \
  return ty;                                                            \
})                                                                      \
template<class dummy=void>                                              \
double ATOMIC_NAME (const double *tx) {                                 \
  double ty[1];                                                         \
  ATOMIC_DOUBLE;                                                        \
  return ty[0];                                                         \
}                                                                       \
template<class dummy=void>                                              \
TMBad::ad_aug ATOMIC_NAME (const TMBad::ad_aug *tx) {                   \
  CppAD::vector<TMBad::ad_aug> tx_(INPUT_SIZE);                         \
  for (size_t i=0; i<INPUT_SIZE; i++) tx_[i]=tx[i];                     \
  return ATOMIC_NAME(tx_)[0];                                           \
}                                                                       \
template<class dummy=void>                                              \
struct ATOMIC_NAME ## Op : TMBad::global::Operator<INPUT_SIZE, 1> {     \
  ATOMIC_NAME ## Op () {}                                               \
  const char* op_name() { return #ATOMIC_NAME; }                        \
  void forward(TMBad::ForwardArgs<TMBad::Scalar> _args_) {              \
    TMBad::Scalar tx[INPUT_SIZE];                                       \
    TMBad::Scalar ty[1]         ;                                       \
    for (size_t i=0; i<INPUT_SIZE; i++) tx[i] = _args_.x(i);            \
    ATOMIC_DOUBLE;                                                      \
    for (size_t i=0; i<1; i++) _args_.y(i) = ty[i];                     \
  }                                                                     \
  static const bool add_forward_replay_copy = true;                     \
  template<class Type> void no_W_set_but_not_used(Type *p) { }          \
  template<class Type> void reverse(TMBad::ReverseArgs<Type> _args_) {  \
    Type tx[INPUT_SIZE];                                                \
    Type ty[1]         ;                                                \
    Type px[INPUT_SIZE];                                                \
    Type py[1]         ;                                                \
    no_W_set_but_not_used(tx);                                          \
    no_W_set_but_not_used(ty);                                          \
    no_W_set_but_not_used(py);                                          \
    for (size_t i=0; i<INPUT_SIZE; i++) tx[i] = _args_.x(i);            \
    for (size_t i=0; i<1         ; i++) ty[i] = _args_.y(i);            \
    for (size_t i=0; i<1         ; i++) py[i] = _args_.dy(i);           \
    ATOMIC_REVERSE;                                                     \
    for (size_t i=0; i<INPUT_SIZE; i++) _args_.dx(i) += px[i];          \
  }                                                                     \
  template<class Type>                                                  \
  void forward                                                          \
  (TMBad::ForwardArgs<Type> &args) { TMBAD_ASSERT(false); }             \
  void reverse                                                          \
  (TMBad::ReverseArgs<TMBad::Writer> &args) { TMBAD_ASSERT(false); }    \
};                                                                      \
template<class dummy>                                                   \
CppAD::vector<TMBad::ad_aug> ATOMIC_NAME                                \
(const CppAD::vector<TMBad::ad_aug> &tx) CSKIP_ATOMIC({                 \
  TMBad::Index m = 1;                                                   \
  typedef ATOMIC_NAME ## Op <> OP;                                      \
  bool all_constant = true;                                             \
  for (size_t i = 0; i<tx.size(); i++)                                  \
    all_constant &= tx[i].constant();                                   \
  CppAD::vector<TMBad::ad_aug> ty(m);                                   \
  if (all_constant) {                                                   \
    CppAD::vector<double> xd(tx.size());                                \
    for (size_t i=0; i<xd.size(); i++) xd[i] = tx[i].Value();           \
    CppAD::vector<double> yd = ATOMIC_NAME(xd);                         \
    for (size_t i=0; i<yd.size(); i++) ty[i] = yd[i];                   \
  } else {                                                              \
    TMBad::OperatorPure*                                                \
      pOp = TMBad::get_glob()->getOperator<OP>();                       \
    std::vector<TMBad::ad_plain>                                        \
      x(&tx[0], &tx[0] + tx.size());                                    \
    std::vector<TMBad::ad_plain>                                        \
      y = TMBad::get_glob()->add_to_stack<OP>(pOp, x);                  \
    for (size_t i=0; i<y.size(); i++) ty[i] = y[i];                     \
  }                                                                     \
  return ty;                                                            \
})                                                                      \
template<class dummy=void>                                              \
void ATOMIC_NAME (const CppAD::vector<TMBad::ad_aug> &tx,               \
                  CppAD::vector<TMBad::ad_aug> &ty) {                   \
  ty = ATOMIC_NAME(tx);                                                 \
}                                                                       \
IF_TMB_PRECOMPILE_ATOMICS(                                              \
template                                                                \
CppAD::vector<double>                                                   \
ATOMIC_NAME<> (const CppAD::vector<double>& tx);                        \
template                                                                \
CppAD::vector<TMBad::ad_aug>                                            \
ATOMIC_NAME<>(const CppAD::vector<TMBad::ad_aug>& tx);                  \
)
// Helper to forward declare atomic
#define TMB_ATOMIC_VECTOR_FUNCTION_DECLARE(ATOMIC_NAME) \
template<class dummy=void>                              \
CppAD::vector<TMBad::ad_aug>                            \
ATOMIC_NAME (const CppAD::vector<TMBad::ad_aug> &x);    \
template<class dummy=void>                              \
CppAD::vector<double>                                   \
ATOMIC_NAME (const CppAD::vector<double> &tx);
/** \brief Construct atomic vector function based on known derivatives

    This macro is used internally to define most atomic functions, see
    for example `atomic::matmul`. It works portably for both CppAD and
    TMBad.

    The macro is actually composed of two other lower-level
    macros. One *declares* and another *defines* the atomic function.
    If two atomic functions depend on one another, the lower-level
    macros may be useful - see e.g the `atomic::fft`.

    \ingroup macros */
#define TMB_ATOMIC_VECTOR_FUNCTION(             \
  ATOMIC_NAME, OUTPUT_DIM,                      \
  ATOMIC_DOUBLE,                                \
  ATOMIC_REVERSE                                \
)                                               \
TMB_ATOMIC_VECTOR_FUNCTION_DECLARE(ATOMIC_NAME) \
TMB_ATOMIC_VECTOR_FUNCTION_DEFINE(              \
  ATOMIC_NAME, OUTPUT_DIM,                      \
  ATOMIC_DOUBLE,                                \
  ATOMIC_REVERSE                                \
)
