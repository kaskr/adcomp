#ifdef TMBAD_FRAMEWORK
/** \brief Sparse and dense versions of atomic Newton solver and Laplace approximation */
namespace newton {
// FIXME: R macro
#ifdef eval
#undef eval
#endif

template <class Type>
struct vector : Eigen::Array<Type, Eigen::Dynamic, 1>
{
  typedef Type value_type;
  typedef Eigen::Array<Type, Eigen::Dynamic, 1> Base;

  vector(void) : Base() {}

  template<class Derived>
  vector(const Eigen::ArrayBase<Derived> &x) : Base(x) {}
  template<class Derived>
  vector(const Eigen::MatrixBase<Derived> &x) : Base(x) {}
  vector(size_t n) : Base(n) {}
  // std::vector
  operator std::vector<Type>() const {
    return std::vector<Type> (Base::data(),
                              Base::data() + Base::size());
  }
  vector(const std::vector<Type> &x) :
    Base(Eigen::Map<const Base> (x.data(), x.size())) { }
};

template <class Type>
struct matrix : Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>
{
  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Base;
  matrix(void) : Base() {}
  template<class Derived>
  matrix(const Eigen::ArrayBase<Derived> &x) : Base(x) {}
  template<class Derived>
  matrix(const Eigen::MatrixBase<Derived> &x) : Base(x) {}
  vector<Type> vec() {
    Base a(*this);
    a.resize(a.size(), 1);
    return a;
  }
};

/** \brief Operator (H, x) -> solve(H, x) */
template <class Hessian_Type>
struct HessianSolveVector : TMBad::global::DynamicOperator< -1, -1 > {
  static const bool have_input_size_output_size = true;
  static const bool add_forward_replay_copy = true;
  /** \warning Pointer */
  Hessian_Type* hessian;
  size_t nnz, x_rows, x_cols; // Dim(x)
  HessianSolveVector(Hessian_Type* hessian, size_t x_cols = 1) :
    hessian ( hessian ),
    nnz     ( hessian->Range() ),
    x_rows  ( hessian->n ),
    x_cols  ( x_cols ) {}
  TMBad::Index input_size() const {
    return nnz + x_rows * x_cols;
  }
  TMBad::Index output_size() const {
    return x_rows * x_cols;
  }
  vector<TMBad::Scalar> solve(const vector<TMBad::Scalar> &h,
                              const vector<TMBad::Scalar> &x) {
    typename Hessian_Type::template MatrixResult<TMBad::Scalar>::type
      H = hessian -> as_matrix(h);
    hessian -> llt_factorize(H); // Assuming analyzePattern(H) has been called once
    matrix<TMBad::Scalar> xm = x.matrix();
    xm.resize(x_rows, x_cols);
    vector<TMBad::Scalar> y = hessian -> llt_solve(H, xm).vec();
    return y;
  }
  vector<TMBad::Replay> solve(const vector<TMBad::Replay> &h,
                              const vector<TMBad::Replay> &x) {
    std::vector<TMBad::ad_plain> hx;
    hx.insert(hx.end(), h.data(), h.data() + h.size());
    hx.insert(hx.end(), x.data(), x.data() + x.size());
    TMBad::global::Complete<HessianSolveVector> Op(*this);
    std::vector<TMBad::ad_plain> ans = Op(hx);
    std::vector<TMBad::ad_aug> ans2(ans.begin(), ans.end());
    return ans2;
  }
  void forward(TMBad::ForwardArgs<TMBad::Scalar> &args) {
    size_t   n = output_size();
    vector<TMBad::Scalar>
      h = args.x_segment(0, nnz),
      x = args.x_segment(nnz, n);
    args.y_segment(0, n) = solve(h, x);
  }
  template <class T>
  void reverse(TMBad::ReverseArgs<T> &args) {
    size_t n = output_size();
    vector<T>
      h  = args. x_segment(0, nnz),
      y  = args. y_segment(0, n),
      dy = args.dy_segment(0, n);
    vector<T> y2 = solve(h, dy);
    for (size_t j=0; j < x_cols; j++) {
      vector<T> y_j  = y .segment(j * x_rows, x_rows);
      vector<T> y2_j = y2.segment(j * x_rows, x_rows);
      vector<T> y2y_j = hessian -> crossprod(y2_j, y_j);
      args.dx_segment(0, nnz) -= y2y_j;
      args.dx_segment(nnz + j * x_rows, x_rows) += y2_j;
    }
  }
  template<class T>
  void forward(TMBad::ForwardArgs<T> &args) { ASSERT(false); }
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) { ASSERT(false); }
  const char* op_name() { return "JSolve"; }
};

/* =================== TODO ===================
   - Optimize RefOp
   - Cleanup get_segment mess
   - Optimize sparsematrix CTOR: Only do 'setFromTriplets' once !
   - Graph permutation: outer - inner trick.
   - SparseMatrix: Diagonal must always be part of pattern !
   - Configurable factorization type - add as template parameter to hessian_types
   - Option for saddlepoint solver ?
   - Move entire code to TMBad ?
   - Interface considerations:
   ---> In addition to solution one wants the hessian evaluated in the optimum. And function value.
   ---> However, the type of the hessian depends on sparse/dense flag.
   - 'hessian_type_dense' rename 'jacobian_dense_t' and similar
   - Generalize HessianSolveVector to handle non-symmetric jacobian
   ============================================
*/

/** \brief Methods specific for a dense hessian */
template<class Factorization=Eigen::LLT<Eigen::Matrix<double, -1, -1> > >
struct jacobian_dense_t : TMBad::ADFun<> {
  typedef TMBad::ADFun<> Base;
  template<class T>
  struct MatrixResult {
    typedef matrix<T> type;
  };
  size_t n;
  Factorization llt;
  jacobian_dense_t() {}
  // FIXME: Want const &G
  // -->   JacFun, var2op, get_keep_var  -->  const
  jacobian_dense_t(TMBad::ADFun<> &H, size_t n) :
    n(n) {
    Base::operator= ( H );
  }
  jacobian_dense_t(TMBad::ADFun<> &F, TMBad::ADFun<> &G, size_t n) :
    n(n) {
    std::vector<bool> keep_x(n, true); // inner
    keep_x.resize(G.Domain(), false);  // outer
    std::vector<bool> keep_y(n, true); // inner
    Base::operator= ( G.JacFun(keep_x, keep_y) );
  }
  template<class V>
  matrix<typename V::value_type> as_matrix(const V &Hx) {
    typedef typename V::value_type T;
    return Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > (Hx.data(), n, n);
  }
  template<class T>
  std::vector<T> eval(const std::vector<T> &x) {
    return Base::operator()(x);
  }
  template<class T>
  matrix<T> operator()(const std::vector<T> &x) {
    return as_matrix(Base::operator()(x));
  }
  /** \brief Sparsity restricted cross product
      In dense case this is the usual cross product
  */
  template<class T>
  vector<T> crossprod(const vector<T> &y2, const vector<T> &y) {
    matrix<T> ans = ( y2.matrix() * y.matrix().transpose() ).array();
    return ans.vec();
  }
  // Sparse.factorize() == Dense.compute()
  void llt_factorize(const matrix<TMBad::Scalar> &h) {
    llt.compute(h);
  }
  Eigen::ComputationInfo llt_info() {
    return llt.info();
  }
  /** \note Optional: This method allows the assumption that a prior
      call to `llt_factorize` has been performed for the same H */
  matrix<TMBad::Scalar> llt_solve(const matrix<TMBad::Scalar> &H,
                                  const matrix<TMBad::Scalar> &x) {
    return llt.solve(x);
  }
  template<class T>
  vector<T> solve(const vector<T> &h,
                  const vector<T> &x) {
    return HessianSolveVector<jacobian_dense_t>(this).solve(h, x);
  }
};
/** \brief Methods specific for a sparse hessian */
template<class Factorization=Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > >
struct jacobian_sparse_t : TMBad::Sparse<TMBad::ADFun<> > {
  typedef TMBad::Sparse<TMBad::ADFun<> > Base;
  template<class T>
  struct MatrixResult {
    typedef Eigen::SparseMatrix<T> type;
  };
  size_t n;
  // FIXME: llt.analyzePattern(H) on construct !!!
  Factorization llt;
  jacobian_sparse_t& operator=(const jacobian_sparse_t &other) {
    Base::operator=(other);
    n = other.n;
    init_llt();
    return *this;
  }
  jacobian_sparse_t (const jacobian_sparse_t &other) : Base(other), n(other.n) {
    init_llt();
  }
  jacobian_sparse_t() {}
  void init_llt() {
    // Analyze pattern
    std::vector<TMBad::Scalar> dummy(this->Range(), 0);
    Eigen::SparseMatrix<TMBad::Scalar> H_dummy = as_matrix(dummy);
    llt.analyzePattern(H_dummy);
  }
  // FIXME: G const !!!
  jacobian_sparse_t(TMBad::Sparse<TMBad::ADFun<> > &H, size_t n) :
    n(n) {
    Base::operator= ( H );
    init_llt();
  }
  jacobian_sparse_t(TMBad::ADFun<> &F, TMBad::ADFun<> &G, size_t n) :
    n(n) {
    std::vector<bool> keep_x(n, true); // inner
    keep_x.resize(G.Domain(), false);  // outer
    std::vector<bool> keep_y(n, true); // inner
    Base::operator= (G.SpJacFun(keep_x, keep_y));
    init_llt();
  }
  template<class V>
  Eigen::SparseMatrix<typename V::value_type> as_matrix(const V &Hx) {
    typedef typename V::value_type T;
    typedef Eigen::Triplet<T> T3;
    std::vector<T3> tripletList;
    size_t K = Hx.size();
    for(size_t k=0; k<K; k++) {
      tripletList.push_back( T3( Base::i[k],
                                 Base::j[k],
                                 Hx[k] ) );
    }
    Eigen::SparseMatrix<T> mat(n, n);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
  }
  template<class T>
  std::vector<T> eval(const std::vector<T> &x) {
    return Base::operator()(x);
  }
  template<class T>
  Eigen::SparseMatrix<T> operator()(const std::vector<T> &x) {
    std::vector<T> Hx = Base::operator()(x);
    return as_matrix(Hx);
  }
  /* Sparsity restricted cross product */
  template<class T>
  vector<T> crossprod(const vector<T> &y2, const vector<T> &y) {
    size_t nnz = Base::Range();
    vector<T> ans(nnz);
    for (size_t k=0; k<nnz; k++) {
      ans[k] = y2[ Base::i [k] ] * y[ Base::j [k] ];
    }
    return ans;
  }
  // Sparse.factorize() == Dense.compute()
  void llt_factorize(const Eigen::SparseMatrix<TMBad::Scalar> &h) {
    llt.factorize(h);
  }
  Eigen::ComputationInfo llt_info() {
    return llt.info();
  }
  /** \note Optional: This method allows the assumption that a prior
      call to `llt_factorize` has been performed for the same H */
  matrix<TMBad::Scalar> llt_solve(const Eigen::SparseMatrix<TMBad::Scalar> &H,
                                  const matrix<TMBad::Scalar> &x) {
    return llt.solve(x);
  }
  template<class T>
  vector<T> solve(const vector<T> &h,
                  const vector<T> &x) {
    return HessianSolveVector<jacobian_sparse_t>(this).solve(h, x);
  }
};

/** \brief Tag operator */
struct TagOp : TMBad::global::Operator<1> {
  static const bool have_eval = true;
  static const bool add_forward_replay_copy = true;
  template<class Type> Type eval(Type x0) {
    return x0 ;
  }
  template<class Type> void reverse(TMBad::ReverseArgs<Type> &args) {
    args.dx(0) += args.dy(0);
  }
  const char* op_name() { return "TagOp"; }
};
/** \brief Mark a variable that links to 'many' random effect */
TMBad::ad_plain Tag(const TMBad::ad_plain &x) {
  return TMBad::get_glob()->add_to_stack<TagOp>(x);
}

/** \brief Methods specific for a sparse plus low rank hessian */
struct jacobian_sparse_plus_lowrank_t {
  // The three tapes
  jacobian_sparse_t<> H;
  TMBad::ADFun<> G;
  jacobian_dense_t<> H0;
  // ADFun methods that should apply to each of the three tapes
  void optimize() {
    H.optimize();
    G.optimize();
    H0.optimize();
  }
  void DomainVecSet(const std::vector<TMBad::Scalar> &x) {
    H.DomainVecSet(x);
    G.DomainVecSet(x);
    H0.DomainVecSet(x);
  }
  void SwapInner() {
    H.SwapInner();
    G.SwapInner();
    H0.SwapInner();
  }
  void SwapOuter() {
    H.SwapOuter();
    G.SwapOuter();
    H0.SwapOuter();
  }
  void print(TMBad::print_config cfg) {
    H.print(cfg);
    G.print(cfg);
    H0.print(cfg);
  }
  // Return type to represent the matrix
  template<class T>
  struct sparse_plus_lowrank {
    Eigen::SparseMatrix<T> H;
    matrix<T> G;
    matrix<T> H0;
    Eigen::Diagonal<Eigen::SparseMatrix<T> > diagonal() {
      return H.diagonal();
    }
  };
  template<class T>
  struct MatrixResult {
    typedef sparse_plus_lowrank<T> type;
  };
  size_t n;
  jacobian_sparse_plus_lowrank_t() {}
  jacobian_sparse_plus_lowrank_t(TMBad::ADFun<> &F,
                                 TMBad::ADFun<> &G_,
                                 size_t n) : n(n) {
    TMBad::Decomp2<TMBad::ADFun<TMBad::ad_aug> >
      F2 = F.decompose("TagOp");
    size_t k = F2.first.Range();
    std::vector<bool> keep_rc(n, true); // inner
    keep_rc.resize(F.Domain(), false);  // outer
    TMBad::Decomp3<TMBad::ADFun<TMBad::ad_aug> >
      F3 = F2.HesFun(keep_rc, true, false, false);
    H = jacobian_sparse_t<>(F3.first, n);
    G = F3.second;
    H0 = jacobian_dense_t<>(F3.third, k);
  }
  // unserialize
  template<class V>
  sparse_plus_lowrank<typename V::value_type> as_matrix(const V &Hx) {
    typedef typename V::value_type T;
    const T* start = Hx.data();
    std::vector<T> v1(start, start + H.Range());
    start += H.Range();
    std::vector<T> v2(start, start + G.Range());
    start += G.Range();
    std::vector<T> v3(start, start + H0.Range());
    sparse_plus_lowrank<T> ans;
    ans.H = H.as_matrix(v1);
    ans.G = vector<T>(v2);
    ans.G.resize(n, v2.size() / n);
    ans.H0 = H0.as_matrix(v3);
    return ans;
  }
  template<class T>
  std::vector<T> eval(const std::vector<T> &x) {
    std::vector<T> ans = H.eval(x);
    std::vector<T> ans2 = G(x);
    std::vector<T> ans3 = H0.eval(x);
    ans.insert(ans.end(), ans2.begin(), ans2.end());
    ans.insert(ans.end(), ans3.begin(), ans3.end());
    return ans;
  }
  template<class T>
  sparse_plus_lowrank<T> operator()(const std::vector<T> &x) {
    return as_matrix(eval(x));
  }
  void llt_factorize(const sparse_plus_lowrank<TMBad::Scalar> &h) {
    H.llt_factorize(h.H);
  }
  // FIXME: Diagonal increments should perhaps be applied to both H and H0.
  Eigen::ComputationInfo llt_info() {
    // Note: As long as diagonal increments are only applied to H this
    // is the relevant info:
    return H.llt_info();
  }
  /** \note Optional: This method allows the assumption that a prior
      call to `llt_factorize` has been performed for the same H */
  matrix<TMBad::Scalar> llt_solve(const sparse_plus_lowrank<TMBad::Scalar> &h,
                                  const matrix<TMBad::Scalar> &x) {
    matrix<TMBad::Scalar> W = H.llt_solve(h.H, h.G); // n x k
    matrix<TMBad::Scalar> M = h.H0.inverse() + h.G.transpose() * W;
    matrix<TMBad::Scalar> y1 = H.llt_solve(h.H, x);
    matrix<TMBad::Scalar> y2 = W * M.ldlt().solve(W.transpose() * x);
    return y1 - y2;
  }
  template<class T>
  vector<T> solve(const vector<T> &h,
                  const vector<T> &x) {
    std::cout << "================ ooooooooooooooooooops !!!\n";
    return x;
  }
  vector<TMBad::Scalar> solve(const vector<TMBad::Scalar> &h,
                              const vector<TMBad::Scalar> &x) {
    sparse_plus_lowrank<TMBad::Scalar> H = as_matrix(h);
    llt_factorize(H);
    return llt_solve(H, x.matrix()).array();
  }
};

/** \brief Newton configuration parameters */
struct newton_config {
  /** \brief Max number of iterations */
  int maxit;
  /** \brief Max number of allowed rejections */
  int max_reject;
  /** \brief Print trace info? */
  int trace;
  /** \brief Convergence tolerance of max gradient component */
  double grad_tol;
  /** \brief Convergence tolerance of consequtive function evaluations (not yet used) */
  double step_tol;
  /** \brief Convergence tolerance of consequtive function evaluations (not yet used) */
  double tol10;
  /** \brief Consider initial guess as invalid if the max gradient component is larger than this number */
  double mgcmax;
  /** \brief Initial step size between 0 and 1 */
  double ustep;
  /** \brief Internal parameter controlling ustep updates */
  double power;
  /** \brief Internal parameter controlling ustep updates */
  double u0;
  /** \brief Use *sparse* as opposed to dense hessian ?
      \details
      Using `sparse=true` for problem that is actually dense have been
      observed to results in a slowdown factor of approximately 3. In
      addition, the dense factorization can be accelerated using the
      `EIGEN_USE_MKL_ALL` preprocessor flag.  On the other hand, using
      `sparse=false` (dense) for a problem that is actually sparse can
      result in much bigger slowdowns.
  */
  bool sparse;
  /** \brief Detect an additional low rank contribution in sparse case? */
  bool lowrank;
  bool decompose;
  /** \brief Behaviour on convergence failure: Report nan-solution ? */
  bool on_failure_return_nan;
  /** \brief Behaviour on convergence failure: Throw warning ?*/
  bool on_failure_give_warning;
  void set_defaults(SEXP x = R_NilValue) {
#define SET_DEFAULT(name, value) set_from_real(x, name, #name, value)
    SET_DEFAULT(maxit, 1000);
    SET_DEFAULT(max_reject, 10);
    SET_DEFAULT(trace, 0);
    SET_DEFAULT(grad_tol, 1e-8);
    SET_DEFAULT(step_tol, 1e-8);
    SET_DEFAULT(tol10, 0.001);
    SET_DEFAULT(mgcmax, 1e+60);
    SET_DEFAULT(ustep, 1);
    SET_DEFAULT(power, 0.5);
    SET_DEFAULT(u0, 1e-04);
    SET_DEFAULT(sparse, false);
    SET_DEFAULT(lowrank, false);
    SET_DEFAULT(decompose, true);
    SET_DEFAULT(on_failure_return_nan, true);
    SET_DEFAULT(on_failure_give_warning, true);
#undef SET_DEFAULT
  }
  newton_config() {
    set_defaults();
  }
  newton_config(SEXP x) {
    set_defaults(x);
  }
  template<class T>
  void set_from_real(SEXP x, T &target, const char* name, double value) {
    SEXP y = getListElement(x, name);
    target = (T) ( y != R_NilValue ? REAL(y)[0] : value );
  }
};
template <class Type>
struct newton_config_t : newton_config {
  newton_config_t() : newton_config() {}
  newton_config_t(SEXP x) : newton_config(x) {}
};

/* Generalized newton solver similar to R function TMB:::newton */
template<class Functor, class Type, class Hessian_Type=jacobian_dense_t<> >
struct NewtonOperator : TMBad::global::SharedDynamicOperator {
  static const bool have_input_size_output_size = true;
  static const bool add_forward_replay_copy = true;
  typedef TMBad::Scalar Scalar;
  typedef TMBad::StdWrap<Functor, vector<TMBad::ad_aug> > FunctorExtend;
  TMBad::ADFun<> function, gradient;
  Hessian_Type hessian;
  // Control convergence
  newton_config cfg;
  // Outer parameters
  std::vector<TMBad::ad_aug> par_outer;
  NewtonOperator(Functor &F, vector<Type> start, newton_config cfg)
    : cfg(cfg)
  {
    // Create tape of the user functor and optimize the tape
    function = TMBad::ADFun<> ( FunctorExtend(F), start);
    function.optimize();
    // The tape may contain expressions that do not depend on the
    // inner parameters (RefOp and expressions that only depend on
    // these). Move such expressions to the parent context?
    if (cfg.decompose) {
      function.decompose_refs();
    }
    // The previous operation does not change the domain vector size.
    size_t n_inner = function.Domain();
    ASSERT(n_inner == (size_t) start.size());
    // Turn remaining references to parent contexts into outer
    // parameters. This operation increases function.Domain()
    par_outer = function.resolve_refs();
    // Mark inner parameter subset in full parameter vector
    std::vector<bool> keep_inner(n_inner, true);
    keep_inner.resize(function.Domain(), false);
    // Gradient function
    gradient = function.JacFun(keep_inner);
    gradient.optimize();
    // Hessian
    hessian = Hessian_Type(function, gradient, n_inner);
    hessian.optimize();
  }
  // Helper to swap inner/outer
  void SwapInner() {
    function.SwapInner();
    gradient.SwapInner();
    hessian.SwapInner();
  }
  void SwapOuter() {
    function.SwapOuter();
    gradient.SwapOuter();
    hessian.SwapOuter();
  }
  // Put it self on tape
  vector<TMBad::ad_aug> add_to_tape() {
    TMBad::global::Complete<NewtonOperator> solver(*this);
    std::vector<TMBad::ad_aug> sol = solver(par_outer);
    // Append outer paramaters to solution
    sol.insert(sol.end(), par_outer.begin(), par_outer.end());
    return sol;
  }
  /* From TMB::newton */
  double phi(double u) {
    return 1. / u - 1.;
  }
  double invphi(double x) {
    return 1. / (x + 1.);
  }
  double increase(double u) {
    return cfg.u0 + (1. - cfg.u0) * std::pow(u, cfg.power);
  }
  double decrease(double u) {
    return (u > 1e-10 ?
            1. - increase(1. - u) :
            (1. - cfg.u0) * cfg.power * u);
  }
  vector<Scalar> x_start; // Cached initial guess
  const char* convergence_fail(const char* msg,
                               vector<Scalar> &x) {
    if (cfg.on_failure_give_warning) {
      if (cfg.trace) {
        Rcout << "Newton convergence failure: " << msg << "\n";
      }
      Rf_warning("Newton convergence failure: %s",
                 msg);
    }
    if (cfg.on_failure_return_nan) {
      x.fill(NAN);
    }
    return msg;
  }
  const char* newton_iterate(vector<Scalar> &x) {
    int reject_counter = 0;
    Scalar f_previous = INFINITY;
    const char* msg = NULL;
    if (x_start.size() == x.size()) {
      Scalar f_x_start = function(x_start)[0];
      Scalar f_x = function(x)[0];
      if ( ! std::isfinite(f_x_start) &&
           ! std::isfinite(f_x) ) {
        return
          convergence_fail("Invalid initial guess", x);
      }
      if (function(x_start)[0] < function(x)[0])
        x = x_start;
    }
    for (int i=0; i < cfg.maxit; i++) {
      vector<Scalar> g = gradient(x);
      Scalar mgc = g.abs().maxCoeff();
      if ( ! std::isfinite(mgc) ||
           mgc > cfg.mgcmax) {
        return
          convergence_fail("Inner gradient had non-finite components", x);
      }
      /* FIXME:
        if (any(!is.finite(g)))
            stop("Newton dropout because inner gradient had non-finite components.")
        if (is.finite(mgcmax) && max(abs(g)) > mgcmax)
            stop("Newton dropout because inner gradient too steep.")
        if (max(abs(g)) < grad.tol)
            return(par)
      */
      if (cfg.trace) std::cout << "mgc=" << mgc << " ";
      if (mgc < cfg.grad_tol) {
        x_start = x;
        if (cfg.trace) std::cout << "\n";
        return msg;
      }
      typename Hessian_Type::template MatrixResult<TMBad::Scalar>::type
        H = hessian(std::vector<Scalar>(x));
      if (cfg.trace) std::cout << "ustep=" << cfg.ustep << " ";
      vector<Scalar> diag_cpy = H.diagonal().array();
      while (true) { // FIXME: Infinite loop
        // H := H + phi * I
        H.diagonal().array() = diag_cpy;
        H.diagonal().array() += phi( cfg.ustep );
        // Try to factorize
        hessian.llt_factorize(H);
        if (hessian.llt_info() == 0) break;
        // H not PD ==> Decrease phi
        cfg.ustep = decrease(cfg.ustep);
      }
      // We now have a PD hessian
      // Let's take a newton step and see if it improves...
      vector<Scalar> x_new =
        x - hessian.llt_solve(H, g).array();
      Scalar f = function(x_new)[0];
      if (std::isfinite(f) &&
          f < f_previous + 1e-8) { // Improvement
        // Accept
        cfg.ustep = increase(cfg.ustep);
        f_previous = f;
        x = x_new;
        reject_counter = 0;
      } else { // No improvement
        // Reject
        cfg.ustep = decrease(cfg.ustep);
        reject_counter ++;
        if (reject_counter > cfg.max_reject)
          return
            convergence_fail("Max number of rejections exceeded", x);
      }
      if (cfg.trace) std::cout << "f=" << f << " ";
      if (cfg.trace) std::cout << "\n";
    }
    return
      convergence_fail("Iteration limit exceeded", x);
  }
  TMBad::Index input_size() const {
    return function.DomainOuter();
  }
  TMBad::Index output_size() const {
    return function.DomainInner(); // Inner dimension
  }
  void forward(TMBad::ForwardArgs<Scalar> &args) {
    size_t n = function.DomainOuter();
    std::vector<Scalar>
      x = args.x_segment(0, n);
    // Set *outer* parameters
    SwapOuter(); // swap
    function.DomainVecSet(x);
    gradient.DomainVecSet(x);
    hessian.DomainVecSet(x);
    SwapOuter(); // swap back
    // Run *inner* iterations
    SwapInner(); // swap
    vector<Scalar> sol = function.DomainVec();
    newton_iterate(sol);
    SwapInner(); // swap back
    args.y_segment(0, sol.size()) = sol;
  }
  template<class T>
  void reverse(TMBad::ReverseArgs<T> &args) {
    vector<T>
      w   = args.dy_segment(0, output_size());
    std::vector<T>
      sol = args. y_segment(0, output_size());
    // NOTE: 'hessian' must have full (inner, outer) vector as input
    size_t n = function.DomainOuter();
    std::vector<T>
      x = args.x_segment(0, n);
    std::vector<T> sol_x = sol; sol_x.insert(sol_x.end(), x.begin(), x.end());
    vector<T> hv = hessian.eval(sol_x);
    vector<T> w2 = - hessian.solve(hv, w);
    vector<T> g = gradient.Jacobian(sol_x, w2);
    args.dx_segment(0, n) += g.tail(n);
  }
  template<class T>
  void forward(TMBad::ForwardArgs<T> &args) { ASSERT(false); }
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) { ASSERT(false); }
  const char* op_name() { return "Newton"; }
  void print(TMBad::global::print_config cfg) {
    Rcout << cfg.prefix << "======== function:\n";
    function.print(cfg);
    Rcout << cfg.prefix << "======== gradient:\n";
    gradient.print(cfg);
    Rcout << cfg.prefix << "======== hessian:\n";
    hessian.print(cfg);
  }
};

template<class Type>
Type log_determinant(const matrix<Type> &H) {
  // FIXME: Depending on TMB atomic
  return atomic::logdet(tmbutils::matrix<Type>(H));
}
template<class Type>
Type log_determinant(const Eigen::SparseMatrix<Type> &H) {
  // FIXME: Tape once for 'reasonable' numeric values - then replay
  // (to avoid unpredictable brancing issues)
  Eigen::SimplicialLDLT< Eigen::SparseMatrix<Type> > ldl(H);
  //return ldl.vectorD().log().sum();
  vector<Type> D = ldl.vectorD();
  return D.log().sum();
}

// Interface
template<class Newton>
void Newton_CTOR_Hook(Newton &F, const vector<double> &start) {
  F.sol = start;
  F.newton_iterate(F.sol);
}
template<class Newton>
void Newton_CTOR_Hook(Newton &F, const vector<TMBad::ad_aug> &start) {
  F.sol = F.add_to_tape();
}
template<class Functor, class Type, class Hessian_Type=jacobian_dense_t<> >
struct NewtonSolver : NewtonOperator<Functor, TMBad::ad_aug, Hessian_Type > {
  typedef NewtonOperator<Functor, TMBad::ad_aug, Hessian_Type > Base;
  typedef typename Hessian_Type::template MatrixResult<Type>::type hessian_t;
  vector<Type> sol; // c(sol, par_outer)
  size_t n; // Number of inner parameters
  NewtonSolver(Functor &F, vector<Type> start, newton_config cfg) :
    Base(F, start, cfg), n(start.size()) {
    Newton_CTOR_Hook(*this, start);
  }
  // Get solution
  operator vector<Type>() const {
    return sol.head(n);
  }
  vector<Type> solution() {
    return sol.head(n);
  }
  Type value() {
    return Base::function(std::vector<Type>(sol))[0];
  }
  hessian_t hessian() {
    return Base::hessian(std::vector<Type>(sol));
  }
  Type Laplace() {
    return
      value() +
      .5 * log_determinant( hessian() ) -
      .5 * log(2. * M_PI) * n;
  }
};

template<class Functor, class Type>
NewtonSolver<Functor,
             Type,
             jacobian_sparse_t<> > NewtonSparse(
                                                Functor &F,
                                                Eigen::Array<Type, Eigen::Dynamic, 1> start,
                                                newton_config cfg = newton_config() ) {
  NewtonSolver<Functor, Type, jacobian_sparse_t<> > ans(F, start, cfg);
  return ans;
}

template<class Functor, class Type>
NewtonSolver<Functor,
             Type,
             jacobian_dense_t<> > NewtonDense(
                                                Functor &F,
                                                Eigen::Array<Type, Eigen::Dynamic, 1> start,
                                                newton_config cfg = newton_config() ) {
  NewtonSolver<Functor, Type, jacobian_dense_t<> > ans(F, start, cfg);
  return ans;
}

template<class Functor, class Type>
NewtonSolver<Functor,
             Type,
             jacobian_sparse_plus_lowrank_t> NewtonSparsePlusLowrank(
                                                                     Functor &F,
                                                                     Eigen::Array<Type, Eigen::Dynamic, 1> start,
                                                                     newton_config cfg = newton_config() ) {
  NewtonSolver<Functor, Type, jacobian_sparse_plus_lowrank_t > ans(F, start, cfg);
  return ans;
}

template<class Functor, class Type>
vector<Type> Newton(Functor &F,
                    Eigen::Array<Type, Eigen::Dynamic, 1> start,
                    newton_config cfg = newton_config() ) {
  if (cfg.sparse) {
    if (! cfg.lowrank)
      return NewtonSparse(F, start, cfg);
    else
      return NewtonSparsePlusLowrank(F, start, cfg);
  }
  else
    return NewtonDense(F, start, cfg);
}

/* Laplace */
template<class Functor, class Type>
Type Laplace(Functor &F,
             Eigen::Array<Type, Eigen::Dynamic, 1> &start,
             newton_config cfg = newton_config() ) {
  if (cfg.sparse) {
    auto opt = NewtonSparse(F, start, cfg);
    start = opt.solution();
    return opt.Laplace();
  } else {
    auto opt = NewtonDense(F, start, cfg);
    start = opt.solution();
    return opt.Laplace();
  }
}

// Usage: slice<>(F, random).Laplace(cfg);
template <class ADFun = TMBad::ADFun<> >
struct slice {
  ADFun &F;
  std::vector<TMBad::Index> random;
  slice(ADFun &F,
        //std::vector<TMBad::ad_aug> x,
        std::vector<TMBad::Index> random) :
    F(F), random(random) {}
  typedef TMBad::ad_aug T;
  std::vector<TMBad::ad_aug> x;
  std::vector<T> operator()(const std::vector<T> &x_random) {
    for (size_t i=0; i<random.size(); i++)
      x[random[i]] = x_random[i];
    return F(x);
  }
  ADFun Laplace_(newton_config cfg = newton_config()) {
    ADFun ans;
    std::vector<double> xd = F.DomainVec();
    x = std::vector<T> (xd.begin(), xd.end());
    ans.glob.ad_start();
    TMBad::Independent(x);
    vector<T> start = TMBad::subset(x, random);
    T y = Laplace(*this, start, cfg);
    y.Dependent();
    ans.glob.ad_stop();
    return ans;
  }
};
TMBad::ADFun<> Laplace_(TMBad::ADFun<> &F,
                        const std::vector<TMBad::Index> &random,
                        newton_config cfg = newton_config() ) CSKIP( {
  slice<> S(F, random);
  return S.Laplace_(cfg);
} )

} // End namespace newton
#endif // TMBAD_FRAMEWORK
