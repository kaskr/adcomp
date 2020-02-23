// Demonstrate adaptive solver of TMBad
#include <TMB.hpp>

struct newton_config {
  int maxit, trace;
  double grad_tol, step_tol, tol10, mgcmax, ustep, power, u0;
  newton_config() :
    maxit(1000), trace(0),
    grad_tol(1e-8),
    step_tol(1e-8),
    tol10(0.001),
    mgcmax(1e+60),
    ustep(1),
    power(0.5),
    u0(1e-04)
  {}
};

/* Generalized newton solver similar to R function TMB:::newton */
template<class Functor, class Type>
struct NewtonSolver : TMBad::global::SharedDynamicOperator {
  static const bool have_input_size_output_size = true;
  static const bool add_forward_replay_copy = true;
  typedef TMBad::Scalar Scalar;
  typedef TMBad::StdWrap<Functor, vector<TMBad::ad_aug> > FunctorExtend;
  // FIXME: For now consider the dense hessian case
  TMBad::ADFun<> function, gradient, hessian;
  // Control convergence
  newton_config cfg;
  NewtonSolver(Functor &F, vector<Type> start, newton_config cfg)
    : cfg(cfg)
  {
    function = TMBad::ADFun<> ( FunctorExtend(F), start);
    function.optimize();
    gradient = function.JacFun();
    gradient.optimize();
    hessian = gradient.JacFun();
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
    std::vector<TMBad::ad_aug> x = function.resolve_refs();
    std::vector<TMBad::ad_aug> y = gradient.resolve_refs();
    std::vector<TMBad::ad_aug> z = hessian.resolve_refs();
    std::vector<TMBad::ad_plain> xyz;
    xyz.insert(xyz.end(), x.begin(), x.end());
    xyz.insert(xyz.end(), y.begin(), y.end());
    xyz.insert(xyz.end(), z.begin(), z.end());
    TMBad::global::Complete<NewtonSolver> solver(*this);
    std::vector<TMBad::ad_plain> ans = solver(xyz);
    std::vector<TMBad::ad_aug> ans2(ans.begin(), ans.end());
    return ans2;
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
    Eigen::LLT<Eigen::Matrix<double, -1, -1> > llt;
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
      vector<Scalar> h = hessian(x);
      if (cfg.trace) std::cout << "ustep=" << cfg.ustep << " ";
      while (true) {
        matrix<Scalar> hm = asMatrix(h, x.size(), x.size());
        hm.diagonal().array() += phi( cfg.ustep );
        llt.compute(hm);
        if (cfg.trace) std::cout << "info=" << llt.info() << " ";
        if (llt.info() == 0) break;
        cfg.ustep = decrease(cfg.ustep);
      }
      // We now have a PD hessian
      // Let's take a newton step and see if it improves...
      matrix<Scalar> step = llt.solve( g.matrix() );
      x = x - step.array();
      Scalar f = function(x)[0];
      if (f < f_previous + 1e-8) {
        // Improvement
        cfg.ustep = increase(cfg.ustep);
        f_previous = f;
      } else {
        // Reject
        cfg.ustep = decrease(cfg.ustep);
        x = x + step.array();
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
  void forward(TMBad::ForwardArgs<T> &args) { ASSERT(false); }
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
    vector<T> h = hessian(sol_z);
    matrix<T> hm = asMatrix(h, sol.size(), sol.size());
    vector<T> w2 = - atomic::matinv(hm) * w.matrix();
    std::vector<T> sol_y = sol; sol_y.insert(sol_y.end(), y.begin(), y.end());
    vector<T> g = gradient.Jacobian(sol_y, w2);
    vector<T> g_outer = g.tail(n2);
    for (size_t i=0; i < (size_t) g_outer.size(); i++) args.dx(i) += g_outer[i];
  }
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) { ASSERT(false); }
  const char* op_name() { return "Newton"; }
};

template<class Functor>
vector<TMBad::ad_aug> Newton(Functor &F, vector<TMBad::ad_aug> start,
                             newton_config cfg = newton_config() ) {
  NewtonSolver<Functor, TMBad::ad_aug > NS(F, start, cfg);
  return NS.add_to_tape();
}
template<class Functor>
vector<double> Newton(Functor &F, vector<double> start,
                      newton_config cfg = newton_config() ) {
  NewtonSolver<Functor, TMBad::ad_aug > NS(F, start, cfg);
  vector<double> x = start;
  NS.newton_iterate(x);
  return x;
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
  DATA_INTEGER(trace);
  DATA_MATRIX(m);
  PARAMETER_VECTOR(x);
  Functor<TMBad::ad_aug> F(m, x);
  newton_config cfg; cfg.trace = trace;
  vector<Type> sol = Newton(F, x, cfg);
  return sol.sum();
}
