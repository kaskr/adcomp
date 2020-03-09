#ifdef TMBAD_FRAMEWORK
/** \brief Sparse and dense versions of atomic Newton solver and Laplace approximation */
namespace newton {

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

/** Methods specific for a dense hessian */
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
  jacobian_dense_t(TMBad::ADFun<> &G, size_t n) :
    n(n) {
    std::vector<bool> keep_x(n, true); // inner
    keep_x.resize(G.Domain(), false);  // outer
    std::vector<bool> keep_y(n, true); // inner
    Base::operator= ( G.JacFun(false, keep_x, keep_y) );
    // FIXME: JacFun proj argument ???
    // rep(keep, n)
    std::vector<bool> keep_dep;
    for (size_t i=0; i<n; i++)
      keep_dep.insert(keep_dep.begin(),
                      keep_x.begin(),
                      keep_x.end());
    this->glob.dep_index = TMBad::subset(this->glob.dep_index,
                                         keep_dep);
  }
  template<class T>
  matrix<T> as_matrix(const std::vector<T> &Hx) {
    return asMatrix(vector<T>(Hx), n, n);
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
  jacobian_sparse_t(TMBad::ADFun<> &G, size_t n) :
    n(n) {
    std::vector<bool> keep_x(n, true); // inner
    keep_x.resize(G.Domain(), false);  // outer
    std::vector<bool> keep_y(n, true); // inner
    Base::operator= (G.SpJacFun(keep_x, keep_y));
    init_llt();
  }
  template<class T>
  Eigen::SparseMatrix<T> as_matrix(const std::vector<T> &Hx) {
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
};

template <class Hessian_Type>
struct HessianSolveVector : TMBad::global::DynamicOperator< -1, -1 > {
  static const bool have_input_size_output_size = true;
  static const bool add_forward_replay_copy = true;
  /** \warning Pointer */
  Hessian_Type* hessian;
  HessianSolveVector(Hessian_Type* hessian) : hessian(hessian) {}
  TMBad::Index input_size() const { return hessian -> Range() + hessian -> n; }
  TMBad::Index output_size() const { return hessian -> n; }
  template<class T>
  std::vector<T> get_segment(TMBad::ForwardArgs<T> &args, size_t from, size_t size) {
    std::vector<T> ans(size);
    for (size_t i=0; i<size; i++) ans[i] = args.x(from + i);
    return ans;
  }
  template<class T>
  std::vector<T> get_segment(TMBad::ReverseArgs<T> &args, size_t from, size_t size) {
    std::vector<T> ans(size);
    for (size_t i=0; i<size; i++) ans[i] = args.x(from + i);
    return ans;
  }
  template<class T>
  std::vector<T> get_segment_y(TMBad::ReverseArgs<T> &args, size_t from, size_t size) {
    std::vector<T> ans(size);
    for (size_t i=0; i<size; i++) ans[i] = args.y(from + i);
    return ans;
  }
  vector<TMBad::Scalar> eval(const std::vector<TMBad::Scalar> &h,
                             const vector<TMBad::Scalar> &x) {
    typename Hessian_Type::template MatrixResult<TMBad::Scalar>::type
      H = hessian -> as_matrix(h);
    hessian -> llt_factorize(H); // Assuming analyzePattern(H) has been called once
    vector<TMBad::Scalar> y = hessian -> llt.solve(x.matrix()).array();
    return y;
  }
  vector<TMBad::Replay> eval(const std::vector<TMBad::Replay> &h,
                             const std::vector<TMBad::Replay> &x) {
    std::vector<TMBad::ad_plain> hx;
    hx.insert(hx.end(), h.begin(), h.end());
    hx.insert(hx.end(), x.begin(), x.end());
    TMBad::global::Complete<HessianSolveVector> solve(*this);
    std::vector<TMBad::ad_plain> ans = solve(hx);
    std::vector<TMBad::ad_aug> ans2(ans.begin(), ans.end());
    return ans2;
  }
  void forward(TMBad::ForwardArgs<TMBad::Scalar> &args) {
    // FIXME: args.get_segment ???
    size_t nnz = hessian -> Range();
    size_t n = hessian -> n;
    std::vector<TMBad::Scalar> h = get_segment(args, 0, nnz);
    vector<TMBad::Scalar> x = get_segment(args, nnz, n);
    vector<TMBad::Scalar> y = eval(h, x);
    for (size_t i=0; i<n; i++) args.y(i) = y[i];
  }
  template <class T>
  void reverse(TMBad::ReverseArgs<T> &args) {

    size_t nnz = hessian -> Range();
    size_t n = hessian -> n;
    std::vector<T> h = get_segment(args, 0, nnz);
    vector<T> y(n), dy(n);
    for (size_t i=0; i<n; i++) {
      y[i] = args.y(i);
      dy[i] = args.dy(i);
    }
    vector<T> y2 = eval(h, dy);
    vector<T> y2y = hessian->crossprod(y2, y);
    for (size_t k=0; k<nnz; k++) {
      args.dx(k) -= y2y[k];
    }
    for (size_t k=0; k<n; k++) {
      args.dx(nnz + k) += y2[k];
    }
  }
  template<class T>
  void forward(TMBad::ForwardArgs<T> &args) { ASSERT(false); }
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) { ASSERT(false); }
  const char* op_name() { return "JSolve"; }
};

/** \brief Newton configuration parameters */
struct newton_config {
  /** \brief Max number of iterations */
  int maxit;
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
  bool decompose;
  /** \brief Behaviour on convergence failure: Report nan-solution ? */
  bool on_failure_return_nan;
  /** \brief Behaviour on convergence failure: Throw warning ?*/
  bool on_failure_give_warning;
  void set_defaults(SEXP x = R_NilValue) {
#define SET_DEFAULT(name, value) set_from_real(x, name, #name, value)
    SET_DEFAULT(maxit, 1000);
    SET_DEFAULT(trace, 0);
    SET_DEFAULT(grad_tol, 1e-8);
    SET_DEFAULT(step_tol, 1e-8);
    SET_DEFAULT(tol10, 0.001);
    SET_DEFAULT(mgcmax, 1e+60);
    SET_DEFAULT(ustep, 1);
    SET_DEFAULT(power, 0.5);
    SET_DEFAULT(u0, 1e-04);
    SET_DEFAULT(sparse, false);
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
    function = TMBad::ADFun<> ( FunctorExtend(F), start);
    function.optimize();
    if (cfg.decompose) {
      function.decompose_refs();
    }
    size_t n_inner = function.Domain();
    ASSERT(n_inner == (size_t) start.size());
    par_outer = function.resolve_refs(); // Increases function.Domain()
    // Mark inner parameter subset
    std::vector<bool> keep_inner(n_inner, true);
    keep_inner.resize(function.Domain(), false);
    // =========================
    // FIXME: inv_inner and inv_outer are no longer set for grad and hess !!!
    // =========================
    // Grad
    gradient = function.JacFun(false, keep_inner);
    gradient.glob.dep_index.resize(n_inner); // FIXME: JacFun should project ???
    gradient.optimize();
    // Hessian
    hessian = Hessian_Type(gradient, n_inner);
    hessian.optimize();
    // FIXME: Hack
    std::vector<bool> keep_outer(keep_inner); keep_outer.flip();
    gradient.inner_inv_index = TMBad::subset(gradient.glob.inv_index, keep_inner);
    gradient.outer_inv_index = TMBad::subset(gradient.glob.inv_index, keep_outer);
    hessian.inner_inv_index = TMBad::subset(hessian.glob.inv_index, keep_inner);
    hessian.outer_inv_index = TMBad::subset(hessian.glob.inv_index, keep_outer);
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
      Rf_warning("Newton convergence failure: %s",
                 msg);
    }
    if (cfg.on_failure_return_nan) {
      x = NAN;
    }
    return msg;
  }
  const char* newton_iterate(vector<Scalar> &x) {
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
        if (hessian.llt.info() == 0) break;
        // H not PD ==> Decrease phi
        cfg.ustep = decrease(cfg.ustep);
      }
      // We now have a PD hessian
      // Let's take a newton step and see if it improves...
      vector<Scalar> x_new =
        x - hessian.llt.solve( g.matrix() ).array();
      Scalar f = function(x_new)[0];
      if (std::isfinite(f) &&
          f < f_previous + 1e-8) { // Improvement
        // Accept
        cfg.ustep = increase(cfg.ustep);
        f_previous = f;
        x = x_new;
      } else { // No improvement
        // Reject
        cfg.ustep = decrease(cfg.ustep);
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
  template<class T>
  std::vector<T> get_segment(TMBad::ForwardArgs<T> &args, size_t from, size_t size) {
    std::vector<T> ans(size);
    for (size_t i=0; i<size; i++) ans[i] = args.x(from + i);
    return ans;
  }
  template<class T>
  std::vector<T> get_segment(TMBad::ReverseArgs<T> &args, size_t from, size_t size) {
    std::vector<T> ans(size);
    for (size_t i=0; i<size; i++) ans[i] = args.x(from + i);
    return ans;
  }
  void forward(TMBad::ForwardArgs<Scalar> &args) {
    size_t n = function.DomainOuter();
    std::vector<Scalar> x = get_segment(args, 0, n);
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
    for (size_t i=0; i < (size_t) sol.size(); i++) args.y(i) = sol[i];
  }
  template<class T>
  void reverse(TMBad::ReverseArgs<T> &args) {
    vector<T> w(output_size());
    for (size_t i=0; i < (size_t) w.size(); i++) w[i] = args.dy(i);
    std::vector<T> sol(output_size());
    for (size_t i=0; i<sol.size(); i++) sol[i] = args.y(i);
    // NOTE: 'hessian' must have full (inner, outer) vector as input
    size_t n = function.DomainOuter();
    std::vector<T> x = get_segment(args, 0, n);
    std::vector<T> sol_x = sol; sol_x.insert(sol_x.end(), x.begin(), x.end());
    HessianSolveVector<Hessian_Type> solve(&hessian);
    std::vector<T> hv = hessian.eval(sol_x);
    vector<T> w2 = - solve.eval(hv, w);
    vector<T> g = gradient.Jacobian(sol_x, w2);
    vector<T> g_outer = g.tail(n);
    for (size_t i=0; i < (size_t) g_outer.size(); i++) args.dx(i) += g_outer[i];
  }
  template<class T>
  void forward(TMBad::ForwardArgs<T> &args) { ASSERT(false); }
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) { ASSERT(false); }
  const char* op_name() { return "Newton"; }
  void print(TMBad::global::print_config cfg) {
    Rcout << cfg.prefix << "======== function:\n";
    function.glob.print(cfg);
    Rcout << cfg.prefix << "======== gradient:\n";
    gradient.glob.print(cfg);
    Rcout << cfg.prefix << "======== hessian:\n";
    hessian.glob.print(cfg);
  }
};

template<class Type>
Type log_determinant(const matrix<Type> &H) {
  return atomic::logdet(H);
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
                                                vector<Type> start,
                                                newton_config cfg = newton_config() ) {
  NewtonSolver<Functor, Type, jacobian_sparse_t<> > ans(F, start, cfg);
  return ans;
}

template<class Functor, class Type>
NewtonSolver<Functor,
             Type,
             jacobian_dense_t<> > NewtonDense(
                                                Functor &F,
                                                vector<Type> start,
                                                newton_config cfg = newton_config() ) {
  NewtonSolver<Functor, Type, jacobian_dense_t<> > ans(F, start, cfg);
  return ans;
}

template<class Functor, class Type>
vector<Type> Newton(Functor &F, vector<Type> start,
                    newton_config cfg = newton_config() ) {
  if (cfg.sparse)
    return NewtonSparse(F, start, cfg);
  else
    return NewtonDense(F, start, cfg);
}

/* Laplace */
template<class Functor, class Type>
Type Laplace(Functor &F, vector<Type> &start,
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
