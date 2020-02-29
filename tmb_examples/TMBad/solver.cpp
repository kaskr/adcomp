// Demonstrate adaptive solver of TMBad
#include <TMB.hpp>

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
    Base::operator= ( G.JacFun() );
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
  /* Sparsity restricted cross product */
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
/** Methods specific for a sparse hessian */
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
    Base::operator= (G.SpJacFun());
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

struct newton_config {
  int maxit, trace;
  double grad_tol, step_tol, tol10, mgcmax, ustep, power, u0;
  bool sparse;
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
  NewtonOperator(Functor &F, vector<Type> start, newton_config cfg)
    : cfg(cfg)
  {
    function = TMBad::ADFun<> ( FunctorExtend(F), start);
    function.optimize();
    gradient = function.JacFun();
    gradient.optimize();
    hessian = Hessian_Type(gradient, start.size());
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
    // Find outer parameters
    std::vector<TMBad::ad_aug> f_par = function.resolve_refs();
    std::vector<TMBad::ad_aug> g_par = gradient.resolve_refs();
    std::vector<TMBad::ad_aug> h_par = hessian.resolve_refs();
    // Append to 'full_par'
    std::vector<TMBad::ad_aug> full_par;
    full_par.insert(full_par.end(), f_par.begin(), f_par.end());
    full_par.insert(full_par.end(), g_par.begin(), g_par.end());
    full_par.insert(full_par.end(), h_par.begin(), h_par.end());
    // Solver: input outer par -> output inner par
    TMBad::global::Complete<NewtonOperator> solver(*this);
    std::vector<TMBad::ad_aug> sol = solver(full_par);
    // Append solution to full_par
    full_par.insert(full_par.end(), sol.begin(), sol.end());
    return full_par;
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
  void newton_iterate(vector<Scalar> &x) {
    if (x_start.size() == x.size()) {
      if (function(x_start)[0] < function(x)[0])
        x = x_start;
    }
    Scalar f_previous = INFINITY;
    for (int i=0; i < cfg.maxit; i++) {
      vector<Scalar> g = gradient(x);
      Scalar mgc = g.abs().maxCoeff();
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
        return;
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
      if (f < f_previous + 1e-8) { // Improvement
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
    Rf_warning("Newton convergence failure");
  }
  TMBad::Index input_size() const {
    size_t n1 = function.DomainOuter();
    size_t n2 = gradient.DomainOuter();
    size_t n3 = hessian.DomainOuter();
    return n1 + n2 + n3;
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
    size_t n1 = function.DomainOuter();
    size_t n2 = gradient.DomainOuter();
    size_t n3 = hessian.DomainOuter();
    std::vector<Scalar> x = get_segment(args,     0, n1);
    std::vector<Scalar> y = get_segment(args,    n1, n2);
    std::vector<Scalar> z = get_segment(args, n1+n2, n3);
    // Set *outer* parameters
    SwapOuter(); // swap
    function.DomainVecSet(x);
    gradient.DomainVecSet(y);
    hessian.DomainVecSet(z);
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
    // FIXME: 'hessian' must have full (inner, outer) vector as input
    size_t n1 = function.DomainOuter();
    size_t n2 = gradient.DomainOuter();
    size_t n3 = hessian.DomainOuter();
    std::vector<T> x = get_segment(args,     0, n1);
    std::vector<T> y = get_segment(args,    n1, n2);
    std::vector<T> z = get_segment(args, n1+n2, n3);
    std::vector<T> sol_z = sol; sol_z.insert(sol_z.end(), z.begin(), z.end());
    HessianSolveVector<Hessian_Type> solve(&hessian);
    std::vector<T> hv = hessian.eval(sol_z);
    vector<T> w2 = - solve.eval(hv, w);
    std::vector<T> sol_y = sol; sol_y.insert(sol_y.end(), y.begin(), y.end());
    vector<T> g = gradient.Jacobian(sol_y, w2);
    vector<T> g_outer = g.tail(n2);
    for (size_t i=0; i < (size_t) g_outer.size(); i++) args.dx(i) += g_outer[i];
  }
  template<class T>
  void forward(TMBad::ForwardArgs<T> &args) { ASSERT(false); }
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) { ASSERT(false); }
  const char* op_name() { return "Newton"; }
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
  vector<Type> sol; // c(f_par, g_par, h_par, sol)
  size_t n; // Number of inner parameters
  NewtonSolver(Functor &F, vector<Type> start, newton_config cfg) :
    Base(F, start, cfg), n(start.size()) {
    // sol = start;
    // Base::newton_iterate(sol);
    Newton_CTOR_Hook(*this, start);
  }
  // Get solution
  operator vector<Type>() const {
    return sol.tail(n);
  }
  vector<Type> solution() {
    return sol.tail(n);
  }
  Type value() {
    // Construct c(sol, f_par)
    size_t n1 = Base::function.DomainOuter();
    vector<Type> f_par_sol(n1 + n);
    f_par_sol << sol.tail(n), sol.head(n1);
    return Base::function(std::vector<Type>(f_par_sol))[0];
  }
  hessian_t hessian() {
    // Construct c(sol, h_par)
    size_t n1 = Base::function.DomainOuter();
    size_t n2 = Base::gradient.DomainOuter();
    size_t n3 = Base::hessian.DomainOuter();
    vector<Type> h_par_sol(n3 + n);
    h_par_sol << sol.tail(n), sol.segment(n1+n2,n3);
    return Base::hessian(std::vector<Type>(h_par_sol));
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


/*
  This class holds the function we wish to minimize

  - It holds a mix of data and parameters
  - We don't want to specify which is which

  - Solution: v^2
*/
template<class Type>
struct Functor {
  matrix<Type> m;
  vector<Type> v;
  Functor(const matrix<Type> &m, const vector<Type> &v) : m(m), v(v) {}
  Type operator()(const vector<Type> &x) {
    v = v * v;
    vector<Type> y = x - v;
    vector<Type> my = m * y;
    return (y * my).sum() + exp(x).sum();
  }
};


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(m);
  PARAMETER_VECTOR(x);
  Functor<TMBad::ad_aug> F(m, x);
  DATA_STRUCT(cfg, newton_config_t);
  Type nll = Laplace(F, x, cfg);
  ADREPORT(x);
  return nll;
}
