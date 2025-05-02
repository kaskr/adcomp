#ifdef TMBAD_FRAMEWORK

/* =================== TODO ===================
   - Move entire code to TMBad ?
   - Option for saddlepoint solver ?
   - Generalize HessianSolveVector to handle non-symmetric jacobian
   - log_determinant of sparse matrix avoid branching issue
   ============================================
*/

/** \brief Construct new sparse matrix of the *same pattern* as input
    \param S SparseMatrix with pattern
    \param x (optional) vector containing new nonzeros
*/
template<class NewType, class OldType>
Eigen::SparseMatrix<NewType> pattern(const Eigen::SparseMatrix<OldType> &S,
                                     std::vector<NewType> x = std::vector<NewType>(0) ) {
  if (S.nonZeros() > 0 && x.size() == 0) {
    x.resize(S.nonZeros());
  }
  return Eigen::Map<const Eigen::SparseMatrix<NewType> > (S.rows(),
                                                          S.cols(),
                                                          S.nonZeros(),
                                                          S.outerIndexPtr(),
                                                          S.innerIndexPtr(),
                                                          x.data(),
                                                          S.innerNonZeroPtr());
}

#ifdef TMBAD_SUPERNODAL
#include "supernodal_inverse_subset.hpp"
#define INVERSE_SUBSET_TEMPLATE Eigen::SupernodalInverseSubset
typedef Eigen::Accessible_CholmodSupernodalLLT DEFAULT_SPARSE_FACTORIZATION;
inline double logDeterminant(const DEFAULT_SPARSE_FACTORIZATION &llt) {
  return llt.logDeterminant();
}
#else
#include "simplicial_inverse_subset.hpp"
#define INVERSE_SUBSET_TEMPLATE Eigen::SimplicialInverseSubset
typedef Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > DEFAULT_SPARSE_FACTORIZATION;
inline double logDeterminant(const DEFAULT_SPARSE_FACTORIZATION &llt) {
  return 2. * llt.matrixL().nestedExpression().diagonal().array().log().sum();
}
#endif

#include <memory>
/** \brief Highly flexible atomic `Newton()` solver and `Laplace()` approximation

    ### Supported features

    - Several hessian structures: Sparse, dense and sparse plus lowrank.
    - Unlimited nesting (Newton solver within newton solver within ...)
    - AD to any order at runtime
    - AD Hessian of Laplace approximation (because inverse subset derivatives are implemented).
    - Sparse hessian of either simplicial or supernodal kind (via preprocessor flag `TMBAD_SUPERNODAL`).
    - CHOLMOD 64 bit integer versions for very large problems (from Eigen version 3.4).
    - Saddlepoint approximation
*/
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
  vector<Type> vec() const {
    Base a(*this);
    a.resize(a.size(), 1);
    return a;
  }
};

/** \brief Operator (H, x) -> solve(H, x)

    Helper operator required to differentiate a newton solver. The
    matrix type `Hessian_Type` can be sparse or dense.
*/
template <class Hessian_Type>
struct HessianSolveVector : TMBad::global::DynamicOperator< -1, -1 > {
  static const bool have_input_size_output_size = true;
  static const bool add_forward_replay_copy = true;
  std::shared_ptr<Hessian_Type> hessian;
  size_t nnz, x_rows, x_cols; // Dim(x)
  HessianSolveVector(std::shared_ptr<Hessian_Type> hessian, size_t x_cols = 1) :
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
  void forward(TMBad::ForwardArgs<T> &args) { TMBAD_ASSERT(false); }
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) { TMBAD_ASSERT(false); }
  const char* op_name() { return "JSolve"; }
};

/* ======================================================================== */
/* === Matrix types that can be used by the Newton solver ================= */
/* ======================================================================== */

/* --- Dense Hessian ------------------------------------------------------ */

/** \brief Methods specific for a dense hessian */
template<class Factorization=Eigen::LLT<Eigen::Matrix<double, -1, -1> > >
struct jacobian_dense_t : TMBad::ADFun<> {
  typedef TMBad::ADFun<> Base;
  template<class T>
  struct MatrixResult {
    typedef matrix<T> type;
  };
  size_t n;
  std::shared_ptr<Factorization> llt;
  jacobian_dense_t() {}
  // FIXME: Want const &F, &G, &H
  // -->   JacFun, var2op, get_keep_var  -->  const
  jacobian_dense_t(TMBad::ADFun<> &H, size_t n) :
    n(n), llt(std::make_shared<Factorization>()) {
    Base::operator= ( H );
  }
  jacobian_dense_t(TMBad::ADFun<> &F, TMBad::ADFun<> &G, size_t n) :
    n(n), llt(std::make_shared<Factorization>()) {
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
    llt->compute(h);
  }
  Eigen::ComputationInfo llt_info() {
    return llt->info();
  }
  /** \note Optional: This method allows the assumption that a prior
      call to `llt_factorize` has been performed for the same H */
  matrix<TMBad::Scalar> llt_solve(const matrix<TMBad::Scalar> &H,
                                  const matrix<TMBad::Scalar> &x) {
    return llt->solve(x);
  }
  template<class T>
  vector<T> solve(std::shared_ptr<jacobian_dense_t> ptr,
                  const vector<T> &h,
                  const vector<T> &x) {
    return HessianSolveVector<jacobian_dense_t>(ptr).solve(h, x);
  }
};

/* --- Sparse Hessian ----------------------------------------------------- */

/** \brief Methods specific for a sparse hessian */
template<class Factorization=DEFAULT_SPARSE_FACTORIZATION >
struct jacobian_sparse_t : TMBad::Sparse<TMBad::ADFun<> > {
  typedef TMBad::Sparse<TMBad::ADFun<> > Base;
  template<class T>
  struct MatrixResult {
    typedef Eigen::SparseMatrix<T> type;
  };
  size_t n;
  std::shared_ptr<Factorization> llt;
  jacobian_sparse_t& operator=(const jacobian_sparse_t &other) {
    Base::operator=(other);
    n = other.n;
    init_llt(); // llt.analyzePattern(H)
    return *this;
  }
  jacobian_sparse_t (const jacobian_sparse_t &other) : Base(other), n(other.n) {
    init_llt(); // llt.analyzePattern(H)
  }
  jacobian_sparse_t() {}
  void init_llt() {
    llt = std::make_shared<Factorization>();
    // Analyze pattern
    std::vector<TMBad::Scalar> dummy(this->Range(), 0);
    Eigen::SparseMatrix<TMBad::Scalar> H_dummy = as_matrix(dummy);
    llt->analyzePattern(H_dummy);
  }
  // FIXME: &F, &G, &H const !!!
  jacobian_sparse_t(TMBad::Sparse<TMBad::ADFun<> > &H, size_t n) :
    n(n) {
    Base::operator= ( H );
    init_llt(); // llt.analyzePattern(H)
  }
  jacobian_sparse_t(TMBad::ADFun<> &F, TMBad::ADFun<> &G, size_t n) :
    n(n) {
    std::vector<bool> keep_x(n, true); // inner
    keep_x.resize(G.Domain(), false);  // outer
    std::vector<bool> keep_y(n, true); // inner
    Base::operator= (G.SpJacFun(keep_x, keep_y));
    init_llt();
  }
  // FIXME: Optimize sparsematrix CTOR by only doing 'setFromTriplets' once?
  template<class V>
  Eigen::SparseMatrix<typename V::value_type> as_matrix(const V &Hx) {
    typedef typename V::value_type T;
    typedef Eigen::Triplet<T> T3;
    std::vector<T3> tripletList(n);
    // Diagonal must always be part of pattern
    for(size_t i=0; i<n; i++) {
      tripletList[i] = T3(i, i, 0);
    }
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
    llt->factorize(h);
  }
  Eigen::ComputationInfo llt_info() {
    return llt->info();
  }
  /** \note Optional: This method allows the assumption that a prior
      call to `llt_factorize` has been performed for the same H */
  matrix<TMBad::Scalar> llt_solve(const Eigen::SparseMatrix<TMBad::Scalar> &H,
                                  const matrix<TMBad::Scalar> &x) {
    return llt->solve(x);
  }
  template<class T>
  vector<T> solve(std::shared_ptr<jacobian_sparse_t> ptr,
                  const vector<T> &h,
                  const vector<T> &x) {
    return HessianSolveVector<jacobian_sparse_t>(ptr).solve(h, x);
  }
};

/* --- Sparse Plus Low Rank Hessian --------------------------------------- */

/** \brief Operator to mark intermediate variables on the tape */
template<class dummy = void>
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
/** \brief Mark a variable during taping

    \details The 'sparse plus low rank' Hessian
    (`jacobian_sparse_plus_lowrank_t`) requires the user to manually
    mark variables used to identify the low rank contribution.

    For any intermediate variable `s(x)` we can write the objective
    function `f(x)` as `f(x,s(x))`. The hessian takes the form

    ```
    hessian( f(x, s(x)) ) = jac(s)^T * f''_yy * jac(s) + R(x)
    ```

    where the first term is referred to as a 'lowrank contribution'
    and `R(x)` is the remainder term. The user must choose `s(x)` in a
    way such that the remainder term becomes a sparse matrix. A good
    heuristic is (which you'll see by expanding the remainder term):

    - Choose `s(x)` that links to *many* random effects
    - Choose `s(x)` such that `hessian(s)` is sparse

    Example (incomplete) of use:
    \code
    PARAMETER_VECTOR(x);
    Type S = x.exp().sum(); // hessian(S) is sparse!
    S = Tag(S);             // Mark S for lowrank contribution
    Type logS = log(S);     // hessian(log(S)) is not sparse!
    return logS;
    \endcode
*/
TMBad::ad_plain Tag(const TMBad::ad_plain &x) CSKIP( {
  return TMBad::get_glob()->add_to_stack<TagOp<> >(x);
} )
/** \brief Otherwise ignore marks  */
TMBad::Scalar Tag(const TMBad::Scalar &x) CSKIP( {
  return x;
} )

/** \brief Methods specific for a sparse plus low rank hessian

    Represents a positive definite matrix of the form

    `H + G * H0 * G^T`

    where

    - H is sparse n-by-n (positive definite)
    - G is a dense n-by-k low rank matrix
    - H0 is dense k-by-k (preferably negative definite)

    To detect this structure one must use the `Tag()` function to select
    k intermediate variables which will be used to decompose the
    computational graph.

    Formulas
    ========

    The following formulas are applied internally to perform the solve
    (Woodbury matrix identity) and determinant (Matrix determinant
    lemma).

    ## Solve

    ```
    W = solve(H, G)
    M = solve(H0) + G^T * W
    solve(H + G * H0 * G^T) = solve(H) - W * solve(M) * W^T
    ```

    ## determinant

    ```
    det(H + G * H0 * G^T) = det(H) * det(H0) * det(M)
    ```

    However, we use the substitution `H0M := H0 * M` for numerical
    robustness. These formulas are valid even when `H0 = 0`:

    ## Solve

    ```
    W <- solve(H, G)
    H0M = I + H0 * t(G) * W
    solve(H + G * H0 * G^T) = solve(H) - W * solve(H0M) * H0 * t(W)
    ```

    ## Determinant

    ```
    det(H + G * H0 * G^T) = det(H) * det(H0M)
    ```
*/
template<class dummy=void>
struct jacobian_sparse_plus_lowrank_t {
  // The three tapes
  std::shared_ptr<jacobian_sparse_t<> > H;
  std::shared_ptr<TMBad::ADFun<> > G;
  std::shared_ptr<jacobian_dense_t<> > H0;
  // ADFun methods that should apply to each of the three tapes
  void optimize() {
    H -> optimize();
    G -> optimize();
    H0 -> optimize();
  }
  void DomainVecSet(const std::vector<TMBad::Scalar> &x) {
    H -> DomainVecSet(x);
    G -> DomainVecSet(x);
    H0 -> DomainVecSet(x);
  }
  void SwapInner() {
    H -> SwapInner();
    G -> SwapInner();
    H0 -> SwapInner();
  }
  void SwapOuter() {
    H -> SwapOuter();
    G -> SwapOuter();
    H0 -> SwapOuter();
  }
  void print(TMBad::print_config cfg) {
    H -> print(cfg);
    G -> print(cfg);
    H0 -> print(cfg);
  }
  // Return type to represent the matrix
  template<class T>
  struct sparse_plus_lowrank {
    Eigen::SparseMatrix<T> H;
    matrix<T> G;
    matrix<T> H0;
    // Optional: Store serialized representation of H
    vector<T> Hvec;
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
    H = std::make_shared<jacobian_sparse_t<> >(F3.first, n);
    G = std::make_shared<TMBad::ADFun<> >(F3.second);
    H0 = std::make_shared<jacobian_dense_t<> >(F3.third, k);
  }
  // unserialize
  template<class V>
  sparse_plus_lowrank<typename V::value_type> as_matrix(const V &Hx) {
    typedef typename V::value_type T;
    const T* start = Hx.data();
    std::vector<T> v1(start, start + H -> Range());
    start += H -> Range();
    std::vector<T> v2(start, start + G -> Range());
    start += G -> Range();
    std::vector<T> v3(start, start + H0 -> Range());
    sparse_plus_lowrank<T> ans;
    ans.H = H -> as_matrix(v1);
    ans.Hvec = v1;
    ans.G = vector<T>(v2);
    ans.G.resize(n, v2.size() / n);
    ans.H0 = H0 -> as_matrix(v3);
    return ans;
  }
  template<class T>
  std::vector<T> eval(const std::vector<T> &x) {
    std::vector<T> ans = H -> eval(x);
    std::vector<T> ans2 = (*G)(x);
    std::vector<T> ans3 = H0 -> eval(x);
    ans.insert(ans.end(), ans2.begin(), ans2.end());
    ans.insert(ans.end(), ans3.begin(), ans3.end());
    return ans;
  }
  template<class T>
  sparse_plus_lowrank<T> operator()(const std::vector<T> &x) {
    return as_matrix(eval(x));
  }
  void llt_factorize(const sparse_plus_lowrank<TMBad::Scalar> &h) {
    H -> llt_factorize(h.H);
  }
  // FIXME: Diagonal increments should perhaps be applied to both H and H0.
  Eigen::ComputationInfo llt_info() {
    // Note: As long as diagonal increments are only applied to H this
    // is the relevant info:
    return H -> llt_info();
  }
  /** \note Optional: This method allows the assumption that a prior
      call to `llt_factorize` has been performed for the same H */
  matrix<TMBad::Scalar> llt_solve(const sparse_plus_lowrank<TMBad::Scalar> &h,
                                  const matrix<TMBad::Scalar> &x) {
    matrix<TMBad::Scalar> W = H -> llt_solve(h.H, h.G); // n x k
    matrix<TMBad::Scalar> H0M = h.H0 * h.G.transpose() * W;
    H0M.diagonal().array() += TMBad::Scalar(1.);
    matrix<TMBad::Scalar> y1 = H -> llt_solve(h.H, x);
    matrix<TMBad::Scalar> y2 = W * H0M.ldlt().solve(h.H0 * W.transpose() * x);
    return y1 - y2;
  }
  template<class T>
  vector<T> solve(std::shared_ptr<jacobian_sparse_plus_lowrank_t<> > ptr,
                  const vector<T> &hvec,
                  const vector<T> &xvec) {
    using atomic::matmul;
    using atomic::matinv;
    sparse_plus_lowrank<T> h = as_matrix(hvec);
    vector<T> s =
      HessianSolveVector<jacobian_sparse_t<> >(ptr -> H,
                                               h.G.cols()).
      solve(h.Hvec, h.G.vec());
    tmbutils::matrix<T> W = s.matrix();
    W.resize(n, W.size() / n);
    tmbutils::matrix<T> H0 = h.H0.array();
    tmbutils::matrix<T> Gt = h.G.transpose();
    tmbutils::matrix<T> H0M =
      matmul(H0,
             matmul(Gt,
                    W));
    H0M.diagonal().array() += T(1.);
    vector<T> y1 =
      HessianSolveVector<jacobian_sparse_t<> >(ptr -> H, 1).
      solve(h.Hvec, xvec);
    tmbutils::matrix<T> iH0M = matinv(H0M);
    tmbutils::matrix<T> Wt = W.transpose();
    tmbutils::matrix<T> xmat = xvec.matrix();
    vector<T> y2 =
      matmul(W,
             matmul(iH0M,
                    matmul(H0,
                           matmul(Wt,
                                  xmat)))).array();
    return y1 - y2;
  }
  // Scalar case: A bit faster than the above template
  vector<TMBad::Scalar> solve(const vector<TMBad::Scalar> &h,
                              const vector<TMBad::Scalar> &x) {
    sparse_plus_lowrank<TMBad::Scalar> H = as_matrix(h);
    llt_factorize(H);
    return llt_solve(H, x.matrix()).array();
  }
  // Helper to get determinant: det(H)*det(H0)*det(M)
  template<class T>
  tmbutils::matrix<T> getH0M(std::shared_ptr<jacobian_sparse_plus_lowrank_t<> > ptr,
                             const sparse_plus_lowrank<T> &h) {
    vector<T> s =
      HessianSolveVector<jacobian_sparse_t<> >(ptr -> H,
                                               h.G.cols()).
      solve(h.Hvec, h.G.vec());
    tmbutils::matrix<T> W = s.matrix();
    W.resize(n, W.size() / n);
    tmbutils::matrix<T> H0 = h.H0.array();
    tmbutils::matrix<T> Gt = h.G.transpose();
    tmbutils::matrix<T> H0M = atomic::matmul(H0, atomic::matmul(Gt, W));
    H0M.diagonal().array() += T(1.);
    return H0M;
  }
};

/* ======================================================================== */
/* === Newton solver ====================================================== */
/* ======================================================================== */

/* --- Configuration ------------------------------------------------------ */

/** \brief Newton configuration parameters */
struct newton_config {
  /** \brief Max number of iterations */
  int maxit;
  /** \brief Max number of allowed rejections */
  int max_reject;
  /** \brief `max_reject` exceeded is convergence success **provided** that Hessian is PD? */
  int ok_exit_if_pdhess;
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
  /** \brief Detect and apply 'dead gradients' simplification */
  bool simplify;
  /** \brief Behaviour on convergence failure: Report nan-solution ? */
  bool on_failure_return_nan;
  /** \brief Behaviour on convergence failure: Throw warning ?*/
  bool on_failure_give_warning;
  /** \brief Consider this absolute reduction 'significant' */
  double signif_abs_reduction;
  /** \brief Consider this relative reduction 'significant' */
  double signif_rel_reduction;
  /** \brief Modify Laplace approximation to return saddlepoint approximation?
      \details
      For this to work the inner objective function must return the
      negative log **MGF** rather than **density** */
  bool SPA;
  void set_defaults(SEXP x = R_NilValue) {
#define SET_DEFAULT(name, value) set_from_real(x, name, #name, value)
    SET_DEFAULT(maxit, 1000);
    SET_DEFAULT(max_reject, 10);
    SET_DEFAULT(ok_exit_if_pdhess, 1);
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
    SET_DEFAULT(simplify, true);
    SET_DEFAULT(on_failure_return_nan, true);
    SET_DEFAULT(on_failure_give_warning, true);
    SET_DEFAULT(signif_abs_reduction, 1e-6);
    SET_DEFAULT(signif_rel_reduction, .5);
    SET_DEFAULT(SPA, false);
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

/* --- Newton Operator ---------------------------------------------------- */

/** \brief Generalized newton solver similar to TMB R function 'newton'
    \details This operator represents a Newton solver that can be put
    on the AD tape. Using the notation
    \f[
    f(u,\theta)
    \f]
    for an objective function with 'inner' parameters \f$u\f$ and
    'outer' parameters \f$\theta\f$, the operator is defined by
    \f[
    \theta \rightarrow \hat{u}(\theta) := \text{argmin}_u f(u, \theta)
    \f].

    ### Applied optimizations

    1. If `cfg.decompose=true` (default), all sub-expressions of
    \f$f(u,\theta)\f$ that do *not* depend on \f$u\f$ are moved out of
    the solver, effecively corresponding to a re-parameterization of
    the outer parameters. An example could be that the objective
    function calculates a determinant that only depends on
    \f$\theta\f$. In this situation the result of the determinant
    becomes the new outer parameter rather than \f$\theta\f$.

    2. A sligthly overlapping, but different, optimization is applied
    if `cfg.simplify=true` (default). We call it a 'dead gradient'
    analysis because it detects which \f$\theta_i\f$ do not affect the
    gradient of the objective wrt \f$u\f$. Such \f$\theta_i\f$ cannot
    affect the solution \f$\hat{u}\f$ either and can therefore be
    removed from the list of outer parameters, effectively reducing the
    input dimension of the Newton operator.

    \tparam Functor Class of the objective function.
    \tparam Hessian_Type Class of the hessian structure (sparse, dense or sparse_plus_lowrank).
*/
template<class Functor, class Hessian_Type=jacobian_dense_t<> >
struct NewtonOperator : TMBad::global::DynamicOperator< -1, -1> {
  static const bool have_input_size_output_size = true;
  static const bool add_forward_replay_copy = true;
  typedef TMBad::Scalar Scalar;
  typedef TMBad::StdWrap<Functor, vector<TMBad::ad_aug> > FunctorExtend;
  TMBad::ADFun<> function, gradient;
  std::shared_ptr<Hessian_Type> hessian;
  // Control convergence
  newton_config cfg;
  // Outer parameters
  std::vector<TMBad::ad_aug> par_outer;
  /** \brief Constructor
      \param F Objective function taking `vector<TMBad::ad_aug>` as input and `TMBad::ad_aug` as output.
      \param start Initial guess for the optimizer.
      \param cfg Configuration parameters - see `newton_config`.
  */
  NewtonOperator(Functor &F, vector<TMBad::ad_aug> start, newton_config cfg)
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
    TMBAD_ASSERT(n_inner == (size_t) start.size());
    // Turn remaining references to parent contexts into outer
    // parameters. This operation increases function.Domain()
    par_outer = function.resolve_refs();
    // Mark inner parameter subset in full parameter vector
    std::vector<bool> keep_inner(n_inner, true);
    keep_inner.resize(function.Domain(), false);
    // Gradient function
    gradient = function.JacFun(keep_inner);
    if (cfg.simplify) {
      // Masks
      std::vector<bool> active = gradient.activeDomain();
      for (size_t i=0; i<n_inner; i++) active[i] = true; // just in case ...
      size_t num_inactive = std::count(active.begin(), active.end(), false);
      if (cfg.trace) {
	std::cout << "Dead gradient args to 'simplify': ";
	std::cout << num_inactive << "\n";
      }
      if (num_inactive > 0) {
	function.DomainReduce(active);
	gradient.DomainReduce(active);
	std::vector<bool> active_outer(active.begin() + n_inner, active.end());
	par_outer = TMBad::subset(par_outer, active_outer);
	TMBAD_ASSERT(n_inner == (size_t) function.inner_inv_index.size());
	function.optimize();
      }
    }
    gradient.optimize();
    // Hessian
    hessian = std::make_shared<Hessian_Type>(function, gradient, n_inner);
    hessian -> optimize();
  }
  // Helper to swap inner/outer
  void SwapInner() {
    function.SwapInner();
    gradient.SwapInner();
    hessian -> SwapInner();
  }
  void SwapOuter() {
    function.SwapOuter();
    gradient.SwapOuter();
    hessian -> SwapOuter();
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
    const char* msg = NULL;
    Scalar f_x = function(x)[0];
    if (x_start.size() == x.size()) {
      Scalar f_x_start = function(x_start)[0];
      if ( ! std::isfinite(f_x_start) &&
           ! std::isfinite(f_x) ) {
        return
          convergence_fail("Invalid initial guess", x);
      }
      if (f_x_start < f_x || ! std::isfinite(f_x)) {
        x = x_start;
        f_x = f_x_start;
      }
    }
    Scalar f_previous = f_x;
    if (cfg.trace) std::cout << "f_start=" << f_x << "\n";
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
        H = (*hessian)(std::vector<Scalar>(x));
      vector<Scalar> diag_cpy = H.diagonal().array();
      // Quick ustep reduction based on Hessian diagonal
      Scalar m = diag_cpy.minCoeff();
      if (std::isfinite(m) && m < 0) {
        Scalar ustep_max = invphi(-m);
        cfg.ustep = std::min(cfg.ustep, ustep_max);
      }
      if (cfg.trace) std::cout << "ustep=" << cfg.ustep << " ";
      while (true) { // FIXME: Infinite loop
        // H := H + phi * I
        H.diagonal().array() = diag_cpy;
        H.diagonal().array() += phi( cfg.ustep );
        // Try to factorize
        hessian -> llt_factorize(H);
        if (hessian -> llt_info() == 0) break;
        // H not PD ==> Decrease phi
        cfg.ustep = decrease(cfg.ustep);
      }
      // We now have a PD hessian
      // Let's take a newton step and see if it improves...
      vector<Scalar> x_new =
        x - hessian -> llt_solve(H, g).array();
      Scalar f = function(x_new)[0];
      // Accept/Reject rules
      bool accept = std::isfinite(f);
      if ( ! accept ) {
        reject_counter++;
      } else {
	// Accept if there's a value improvement:
	accept =
	  (f < f_previous);
	// OR...
	if (! accept) {
	  // Accept if relative mgc reduction is substantial without increasing value 'too much'
	  accept =
	    (f - f_previous <= cfg.signif_abs_reduction ) &&
	    vector<Scalar>(gradient(x_new)).abs().maxCoeff() / mgc < 1. - cfg.signif_rel_reduction;
	}
	// Assess the quality of the update
	if ( accept && (f_previous - f > cfg.signif_abs_reduction) ) {
	  reject_counter = 0;
	} else {
	  reject_counter++;
	}
      }
      // Take action depending on accept/reject
      if ( accept ) { // Improvement
        // Accept
        cfg.ustep = increase(cfg.ustep);
        f_previous = f;
        x = x_new;
      } else { // No improvement
        // Reject
        cfg.ustep = decrease(cfg.ustep);
      }
      // Handle long runs without improvement
      if (reject_counter > cfg.max_reject) {
	if (cfg.ok_exit_if_pdhess) {
	  H = (*hessian)(std::vector<Scalar>(x));
	  // Try to factorize
	  hessian -> llt_factorize(H);
	  bool PD = (hessian -> llt_info() == 0);
	  if (cfg.trace) std::cout << "Trying early exit - PD Hess? " << PD << "\n";
	  if (PD) return msg;
	}
	return
	  convergence_fail("Max number of rejections exceeded", x);
      }
      // Tracing info
      if (cfg.trace) std::cout << "f=" << f << " ";
      if (cfg.trace) std::cout << "reject=" << reject_counter << " ";
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
    hessian -> DomainVecSet(x);
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
    vector<T> hv = hessian -> eval(sol_x);
    vector<T> w2 = - hessian -> solve(hessian, hv, w);
    vector<T> g = gradient.Jacobian(sol_x, w2);
    args.dx_segment(0, n) += g.tail(n);
  }
  template<class T>
  void forward(TMBad::ForwardArgs<T> &args) { TMBAD_ASSERT(false); }
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) { TMBAD_ASSERT(false); }
  const char* op_name() { return "Newton"; }
  void print(TMBad::global::print_config cfg) {
    Rcout << cfg.prefix << "======== function:\n";
    function.print(cfg);
    Rcout << cfg.prefix << "======== gradient:\n";
    gradient.print(cfg);
    Rcout << cfg.prefix << "======== hessian:\n";
    hessian -> print(cfg);
  }
};

typedef jacobian_sparse_t< DEFAULT_SPARSE_FACTORIZATION > jacobian_sparse_default_t;

template<class Factorization=DEFAULT_SPARSE_FACTORIZATION >
struct InvSubOperator : TMBad::global::DynamicOperator< -1, -1> {
  static const bool have_input_size_output_size = true;
  static const bool add_forward_replay_copy = true;
  typedef TMBad::Scalar Scalar;
  Eigen::SparseMatrix<Scalar> hessian; // Pattern
  std::shared_ptr< Factorization > llt; // Factorization
  INVERSE_SUBSET_TEMPLATE<double> ihessian;
  //typedef Eigen::AutoDiffScalar<Eigen::Array<double, 1, 1> > ad1;
  typedef atomic::tiny_ad::variable<1,1> ad1;
  INVERSE_SUBSET_TEMPLATE<ad1> D_ihessian;
  InvSubOperator(const Eigen::SparseMatrix<Scalar> &hessian,
                 std::shared_ptr< Factorization > llt) :
    hessian(hessian), llt(llt), ihessian(llt), D_ihessian(llt) { }
  TMBad::Index input_size() const {
    return hessian.nonZeros();
  }
  TMBad::Index output_size() const {
    return hessian.nonZeros();
  }
  void forward(TMBad::ForwardArgs<Scalar> &args) {
    size_t n = input_size();
    std::vector<Scalar>
      x = args.x_segment(0, n);
    Eigen::SparseMatrix<Scalar> h = pattern(hessian, x);
    llt->factorize(h); // Update shared factorization
    h = ihessian(h);
    args.y_segment(0, n) = h.valuePtr();
  }
  void reverse(TMBad::ReverseArgs<Scalar> &args) {
    size_t n = input_size();
    std::vector<Scalar> x  = args. x_segment(0, n);
    std::vector<Scalar> dy = args.dy_segment(0, n);
    Eigen::SparseMatrix<Scalar> dy_mat = pattern(hessian, dy);
    // Symmetry correction (range direction)
    dy_mat.diagonal() *= 2.; dy_mat *= .5;
    std::vector<ad1> x_ (n);
    for (size_t i=0; i<n; i++) {
      x_[i].value = x[i];
      x_[i].deriv[0] = dy_mat.valuePtr()[i];
    }
    Eigen::SparseMatrix<ad1> h_ = pattern(hessian, x_);
    // Inverse subset
    h_ = D_ihessian(h_);
    // Symmetry correction
    h_.diagonal() *= .5; h_ *= 2.;
    // Gather result
    std::vector<Scalar> dx(n);
    for (size_t i=0; i<n; i++) {
      dx[i] = h_.valuePtr()[i].deriv[0];
    }
    args.dx_segment(0, n) += dx;
  }
  template<class T>
  void reverse(TMBad::ReverseArgs<T> &args) {
    Rf_error("Inverse subset: order 2 not yet implemented (try changing config())");
  }
  template<class T>
  void forward(TMBad::ForwardArgs<T> &args) { TMBAD_ASSERT(false); }
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) { TMBAD_ASSERT(false); }
  const char* op_name() { return "InvSub"; }
};
template<class Factorization=DEFAULT_SPARSE_FACTORIZATION >
struct LogDetOperator : TMBad::global::DynamicOperator< -1, -1> {
  static const bool have_input_size_output_size = true;
  static const bool add_forward_replay_copy = true;
  typedef TMBad::Scalar Scalar;
  Eigen::SparseMatrix<Scalar> hessian; // Pattern
  std::shared_ptr< Factorization > llt; // Factorization
  INVERSE_SUBSET_TEMPLATE<double> ihessian;
  LogDetOperator(const Eigen::SparseMatrix<Scalar> &hessian,
                 std::shared_ptr< Factorization > llt) :
    hessian(hessian), llt(llt), ihessian(llt) { }
  TMBad::Index input_size() const {
    return hessian.nonZeros();
  }
  TMBad::Index output_size() const {
    return 1;
  }
  void forward(TMBad::ForwardArgs<Scalar> &args) {
    size_t n = input_size();
    std::vector<Scalar>
      x = args.x_segment(0, n);
    Eigen::SparseMatrix<Scalar> h = pattern(hessian, x);
    llt->factorize(h);
    // Get out if factorization went wrong
    if (llt->info() != 0) {
      args.y(0) = R_NaN;
      return;
    }
    args.y(0) = logDeterminant(*llt);
  }
  void reverse(TMBad::ReverseArgs<Scalar> &args) {
    size_t n = input_size();
    // Get out if factorization went wrong
    if (llt->info() != 0) {
      for (size_t i=0; i<n; i++) args.dx(i) += R_NaN;
      return;
    }
    std::vector<Scalar> x = args.x_segment(0, n);
    Eigen::SparseMatrix<Scalar> ih = pattern(hessian, x);
    // Inverse subset
    ih = ihessian(ih);
    // Symmetry correction
    ih.diagonal() *= .5; ih *= 2.;
    ih *= args.dy(0);
    args.dx_segment(0, n) += ih.valuePtr();
  }
  void reverse(TMBad::ReverseArgs<TMBad::ad_aug> &args) {
    size_t n = input_size();
    TMBad::global::Complete< InvSubOperator<Factorization> > IS(hessian, llt);
    std::vector<TMBad::ad_aug> x = args.x_segment(0, n);
    std::vector<TMBad::ad_aug> y = IS(x);
    Eigen::SparseMatrix<TMBad::ad_aug> ih = pattern(hessian, y);
    // Symmetry correction
    ih.diagonal() *= .5; ih *= 2.;
    ih *= args.dy(0);
    args.dx_segment(0, n) += ih.valuePtr();
  }
  template<class T>
  void forward(TMBad::ForwardArgs<T> &args) { TMBAD_ASSERT(false); }
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) { TMBAD_ASSERT(false); }
  const char* op_name() { return "logDet"; }
};
template<class Type>
Type log_determinant_simple(const Eigen::SparseMatrix<Type> &H) {
  // FIXME: Tape once for 'reasonable' numeric values - then replay
  // (to avoid unpredictable brancing issues)
  Eigen::SimplicialLDLT< Eigen::SparseMatrix<Type> > ldl(H);
  //return ldl.vectorD().log().sum();
  vector<Type> D = ldl.vectorD();
  return D.log().sum();
}
template<class Type>
Type log_determinant(const Eigen::SparseMatrix<Type> &H,
                     std::shared_ptr< jacobian_sparse_default_t > ptr) {
  if (!config.tmbad.atomic_sparse_log_determinant)
    return log_determinant_simple(H);
  const Type* vptr = H.valuePtr();
  size_t n = H.nonZeros();
  std::vector<Type> x(vptr, vptr + n);
  TMBad::global::Complete<LogDetOperator<> > LD(pattern<double>(H), ptr->llt);
  std::vector<Type> y = LD(x);
  return y[0];
}
template<class Type>
Type log_determinant(const Eigen::SparseMatrix<Type> &H) {
  const Type* vptr = H.valuePtr();
  size_t n = H.nonZeros();
  std::vector<Type> x(vptr, vptr + n);
  Eigen::SparseMatrix<double> H_pattern = pattern<double>(H);
  std::shared_ptr< DEFAULT_SPARSE_FACTORIZATION > llt =
    std::make_shared< DEFAULT_SPARSE_FACTORIZATION > (H_pattern);
  TMBad::global::Complete<LogDetOperator<> > LD(H_pattern, llt);
  std::vector<Type> y = LD(x);
  return y[0];
}
inline double log_determinant(const Eigen::SparseMatrix<double> &H) {
  DEFAULT_SPARSE_FACTORIZATION llt(H);
  return logDeterminant(llt);
}
inline double log_determinant(const Eigen::SparseMatrix<double> &H,
                              std::shared_ptr<jacobian_sparse_default_t> ptr) {
  // FIXME: Use ptr llt
  DEFAULT_SPARSE_FACTORIZATION llt(H);
  return logDeterminant(llt);
}

template<class Type, class PTR>
Type log_determinant(const matrix<Type> &H, PTR ptr) {
  // FIXME: Depending on TMB atomic
  return atomic::logdet(tmbutils::matrix<Type>(H));
}
template<class Type>
Type log_determinant(const jacobian_sparse_plus_lowrank_t<>::sparse_plus_lowrank<Type> &H,
                     std::shared_ptr<jacobian_sparse_plus_lowrank_t<> > ptr) {
  matrix<Type> H0M = (ptr -> getH0M(ptr, H)).array();
  return
    log_determinant(H.H, ptr->H) +
    log_determinant(H0M, NULL);
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

/** Helper to simplify user interface of Functor by making `ad_aug` version work for `double` as well

    In principle, the double case should be trivial, as it evaluates
    the functor using AD types that are constants. However, some
    operators might add to tape for constant inputs (by mistake). To
    guard aginst such mistakes, we wrap the double evaluation inside a
    dummy AD context.
*/
template<class Functor, class Type>
struct safe_eval {
  Type operator()(Functor &F, vector<Type> x) {
    return F(x);
  }
};
template<class Functor>
struct safe_eval<Functor, double> {
  double operator()(Functor &F, vector<double> x) {
    /* double case: Wrap evaluation inside an AD context, just in case
       some operations in F adds to tape. */
    TMBad::global dummy;
    dummy.ad_start();
    double ans = asDouble(F(x));
    dummy.ad_stop();
    return ans;
  }
};

template<class Functor, class Type, class Hessian_Type=jacobian_dense_t<> >
struct NewtonSolver : NewtonOperator<Functor, Hessian_Type > {
  typedef NewtonOperator<Functor, Hessian_Type > Base;
  typedef typename Hessian_Type::template MatrixResult<Type>::type hessian_t;
  vector<Type> sol; // c(sol, par_outer)
  size_t n; // Number of inner parameters
  Functor& F;
  NewtonSolver(Functor &F, vector<Type> start, newton_config cfg) :
    Base(F, start, cfg), n(start.size()), F(F) {
    Newton_CTOR_Hook(*this, start);
  }
  // Get solution
  operator vector<Type>() const {
    return sol.head(n);
  }
  tmbutils::vector<Type> solution() {
    return sol.head(n);
  }
  Type value() {
    if (Base::cfg.simplify) {
      return safe_eval<Functor, Type>()(F, solution());
    } else {
      return Base::function(std::vector<Type>(sol))[0];
    }
  }
  hessian_t hessian() {
    return (*(Base::hessian))(std::vector<Type>(sol));
  }
  Type Laplace() {
    double sign = (Base::cfg.SPA ? -1 : 1);
    return
      sign * value() +
      .5 * log_determinant( hessian(),
                            Base::hessian) -
      sign * .5 * log(2. * M_PI) * n;
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
             jacobian_sparse_plus_lowrank_t<> > NewtonSparsePlusLowrank(
                                                                     Functor &F,
                                                                     Eigen::Array<Type, Eigen::Dynamic, 1> start,
                                                                     newton_config cfg = newton_config() ) {
  NewtonSolver<Functor, Type, jacobian_sparse_plus_lowrank_t<> > ans(F, start, cfg);
  return ans;
}

/** \brief Tape a functor and return solution

    Can be used anywhere in a template. Inner and outer parameters are
    automatically detected.

    \param F Function to minimize
    \param start Vector with initial guess
    \param cfg Configuration parameters for solver

    \return Vector with solution
*/
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

/** \brief Tape a functor and return Laplace Approximation

    Can be used anywhere in a template. Inner and outer parameters are
    automatically detected.

    \param F Function to minimize
    \param start Vector with initial guess
    \param cfg Configuration parameters for solver

    \return Scalar with Laplace Approximation

    \note `start` is passed by reference and contains the inner
    problem mode on output
*/
template<class Functor, class Type>
Type Laplace(Functor &F,
             Eigen::Array<Type, Eigen::Dynamic, 1> &start,
             newton_config cfg = newton_config() ) {
  if (cfg.sparse) {
    if (!cfg.lowrank) {
      auto opt = NewtonSparse(F, start, cfg);
      start = opt.solution();
      return opt.Laplace();
    } else {
      auto opt = NewtonSparsePlusLowrank(F, start, cfg);
      start = opt.solution();
      return opt.Laplace();
    }
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
    F(F), random(random) {
    TMBAD_ASSERT2(F.Range() == 1,
            "Laplace approximation is for scalar valued functions");
  }
  typedef TMBad::ad_aug T;
  std::vector<TMBad::ad_aug> x;
  T operator()(const std::vector<T> &x_random) {
    for (size_t i=0; i<random.size(); i++)
      x[random[i]] = x_random[i];
    return F(x)[0];
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
