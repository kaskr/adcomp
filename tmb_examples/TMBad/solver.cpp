// Demonstrate adaptive solver of TMBad
#include <TMB.hpp>

/* Generalized newton solver similar to R function TMB:::newton */
template<class Functor, class Type>
struct NewtonSolver : TMBad::global::DynamicOperator< -1, -1 > {
  static const bool have_input_size_output_size = true;
  static const bool add_forward_replay_copy = true;
  typedef TMBad::Scalar Scalar;
  typedef TMBad::StdWrap<Functor, vector<TMBad::ad_aug> > FunctorExtend;
  // FIXME: For now consider the dense hessian case
  TMBad::ADFun<> function, gradient, hessian;
  // Control convergence
  int maxit;
  Scalar grad_tol;
  NewtonSolver(Functor &F, vector<Type> start) :
    maxit(100), grad_tol(1e-8)
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
  void newton_iterate(std::vector<Scalar> &x) {
    for (int i=0; i<maxit; i++) {
      vector<Scalar> g = gradient(x);
      Scalar mgc = g.abs().maxCoeff();
      if (mgc < grad_tol) return;
      vector<Scalar> h = hessian(x);
      matrix<Scalar> hm = asMatrix(h, x.size(), x.size());
      vector<Scalar> step = hm.inverse() * g.matrix();
      x = vector<Scalar>(vector<Scalar>(x) - step);
    }
    Rf_warning("Convergence failure");
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
    std::vector<Scalar> x = get_segment(args,  0, n1);
    std::vector<Scalar> y = get_segment(args, n1, n2);
    std::vector<Scalar> z = get_segment(args, n2, n3);
    // Set *outer* parameters
    SwapOuter(); // swap
    function.DomainVecSet(x);
    gradient.DomainVecSet(y);
    hessian.DomainVecSet(z);
    SwapOuter(); // swap back
    // Run *inner* iterations
    SwapInner(); // swap
    std::vector<Scalar> sol = function.DomainVec();
    newton_iterate(sol);
    SwapInner(); // swap back
    for (size_t i=0; i<sol.size(); i++) args.y(i) = sol[i];
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
    std::vector<T> x = get_segment(args,  0, n1);
    std::vector<T> y = get_segment(args, n1, n2);
    std::vector<T> z = get_segment(args, n2, n3);
    std::vector<T> sol_z = sol; sol_z.insert(sol_z.end(), z.begin(), z.end());
    vector<T> h = hessian(sol_z);
    matrix<T> hm = asMatrix(h, sol.size(), sol.size());
    vector<T> w2 = - hm.inverse() * w.matrix();
    std::vector<T> sol_y = sol; sol_y.insert(sol_y.end(), y.begin(), y.end());
    vector<T> g = gradient.Jacobian(sol_y, w2);
    vector<T> g_outer = g.tail(n2);
    for (size_t i=0; i < (size_t) g_outer.size(); i++) args.dx(i) += g_outer[i];
  }
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) { ASSERT(false); }
  const char* op_name() { return "Newton"; }
};

template<class Functor>
vector<TMBad::ad_aug> Newton(Functor &F, vector<TMBad::ad_aug> start) {
  NewtonSolver<Functor, TMBad::ad_aug > NS(F, start);
  return NS.add_to_tape();
}
template<class Functor>
vector<double> Newton(Functor &F, vector<double> start) {
  NewtonSolver<Functor, TMBad::ad_aug > NS(F, start);
  std::vector<double> x = start;
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
    return (y * my).sum();
  }
};


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(m);
  PARAMETER_VECTOR(x);
  Functor<TMBad::ad_aug> F(m, x);
  vector<Type> sol = Newton(F, x);
  return sol.sum();
}
