#ifdef TMBAD_FRAMEWORK

/** \file
    \brief Atomic sparse matrix exponentiation
*/

/** \brief Sparse matrix exponentiation

    This namespace contains methods to calculate `x^T*exp(A)` when A
    is a sparse matrix.
*/
namespace sparse_matrix_exponential {

/* Shared data and methods for sparse matrix vector multiply: A^T * x */
template <class Type>
struct SparseATx {
  // Eigens index type
  typedef typename
  Eigen::SparseMatrix<Type>::StorageIndex StorageIndex;
  // Sparse matrix (A) dimension and number of non zeros
  StorageIndex nrow, ncol, nnz;
  // Sparsity pattern (A)
  std::vector<StorageIndex> i; // Length nnz
  std::vector<StorageIndex> p; // Length ncol+1
  SparseATx(const Eigen::SparseMatrix<Type> &x) :
    nrow(x.rows()),
    ncol(x.cols()),
    nnz(x.nonZeros()),
    i(x.innerIndexPtr(), x.innerIndexPtr() + nnz),
    p(x.outerIndexPtr(), x.outerIndexPtr() + ncol + 1) {}
  StorageIndex rows() const { return nrow; }
  StorageIndex cols() const { return ncol; }
  /* Operation: y = A^T * x */
  template<class T>
  void f(const T *A, const T *x, T *y) {
    for (StorageIndex j = 0; j < ncol; j++) {
      y[j] = 0;
      for (StorageIndex k = p[j]; k < p[j+1]; k++) {
        y[j] += A[k] * x[i[k]]; // good access
      }
    }
  }
  /* Adjoint operation */
  template<class T>
  void df(const T *A, const T *x, const T *y,
          T *dA, T *dx, const T *dy) {
    for (StorageIndex j = 0; j < ncol; j++) {
      for (StorageIndex k = p[j]; k < p[j+1]; k++) {
        dA[k] += x[i[k]] * dy[j]; // good access
        dx[i[k]] += A[k] * dy[j]; // bad access
      }
    }
  }
};


template<class Type, bool transpose=false>
struct SpAxOp : TMBad::global::DynamicOperator< -1, -1> {
  typedef typename
  Eigen::SparseMatrix<Type>::StorageIndex StorageIndex;
  std::shared_ptr<SparseATx<Type> > P;
  SpAxOp() {}
  SpAxOp(const Eigen::SparseMatrix<Type> &x) :
    P(std::make_shared<SparseATx<Type> >( x )) { }
  static const bool have_input_size_output_size = true;
  TMBad::Index input_size() const {
    return 2; // Two pointers
  }
  TMBad::Index output_size() const {
    return (*P).cols();
  }
  // FIXME:
  //   static const bool forward_replay_copy = true;
  template<class T>
  void forward(TMBad::ForwardArgs<T> &args) {
    const T* A = args.x_ptr(0);
    const T* x = args.x_ptr(1);
    T* y = args.y_ptr(0);
    (*P).template f<T>(A, x, y);
  }
  template<class T>
  void reverse(TMBad::ReverseArgs<T> &args) {
    const T* A = args.x_ptr(0);
    const T* x = args.x_ptr(1);
    const T* y = args.y_ptr(0);
    T* dA = args.dx_ptr(0);
    T* dx = args.dx_ptr(1);
    const T* dy = args.dy_ptr(0);
    (*P).template df<T>(A, x, y, dA, dx, dy);
  }
  // ---- Dependencies ---- (copied from ad_blas)
  // FIXME: Make general implementation for pointer based case
  void dependencies(TMBad::Args<> &args, TMBad::Dependencies &dep) const {
    dep.add_segment(args.input(0), (*P).nnz  );
    dep.add_segment(args.input(1), (*P).rows() );
  }
  static const bool have_dependencies = true;
  /** \brief This operator **has** implicit dependencies */
  static const bool implicit_dependencies = true;
  /** \brief It is **not* safe to remap the inputs of this operator */
  static const bool allow_remap = false;
  // Not yet implemented
  void forward(TMBad::ForwardArgs<TMBad::Writer> &args) { ASSERT(false); }
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) { ASSERT(false); }
  const char* op_name() { return "SpAxOp"; }
};

// FIXME: Remember to transpose A on construction

/** \brief Shared configuration parameters for `expm_series` and `expm_generator` */
template<class T>
struct config {
  /** \brief Normalize result? (expm_generator only) */
  bool normalize;
  /** \brief Trace retaping? */
  bool trace;
  /** \brief Warn if Nmax exceeded? */
  bool warn;
  /** \brief Use no more than this number of terms in series. */
  int Nmax;
  config() : normalize(false), trace(true), warn(true), Nmax(100) {}
};

/** \brief Matrix exponential multiplied by vector

    This shared atomic operator can be used in cases where
    `x^T*exp(A)` needs to be calculated many times for same sparsity
    pattern. The number of terms (N) in series can be variable and
    retaping is triggered if (and only if) N changes compared to
    previous evaluation.

    \note Multiplication with vector is from the left: `x^T * exp(A)`
*/
template <class T>
struct expm_series {
  typedef vectorize::vector<T> vec;
  /** Number of terms */
  T N;
  /** Vector of sparse matrix non zeros */
  vec A_values;
  /** Operator that multiplies sparse matrix with vector */
  TMBad::global::Complete<SpAxOp<T> > multiply;
  /** ADFun object that holds the entire matexp tape */
  TMBad::ADFun_packed<> F;
  /** Configuration */
  config<T> cfg;
  /** Helper to update generator */
  void update(Eigen::SparseMatrix<T> &A) {
    // FIXME: Assert same pattern
    A_values = vec(A.valuePtr(), A.nonZeros());
  }
  expm_series() {}
  expm_series(Eigen::SparseMatrix<T> &A, T N, config<T> cfg=config<T>()) :
    N(N), A_values(A.valuePtr(), A.nonZeros()), multiply(A), cfg(cfg)
  { }
  /** \brief Evaluate `x^T*exp(A)` */
  vec operator()(vec x) {
    N = min(N, T(cfg.Nmax));
    std::vector<TMBad::ad_segment> args = {A_values, x, vec(N, 1)};
    if (! F.initialized() ) {
      struct Test {
        config<T> cfg;
        TMBad::Scalar Nold;
        bool operator() (const std::vector<TMBad::Scalar*> &x) {
          using TMBad::operator<<;
          TMBad::Scalar N = x[2][0];
          if ( (int) N == cfg.Nmax) {
            if (cfg.warn)
              Rf_warning("expm: N terms reduced to Nmax (%i)", (int) cfg.Nmax);
          }
          bool change = (Nold != N);
          if (cfg.trace && change) {
            Rcout << "Retaping:" << " Nold=" << Nold << " Nnew=" << N << "\n";
            Nold = N;
          }
          return change;
        }
      };
      Test N_changed = {cfg, N.Value()};
      F = TMBad::ADFun_retaping(*this, args, N_changed);
    }
    return F(args);
  }
private:
  friend class TMBad::PackWrap<expm_series>;
  // Packed: (A, x) -> exp(A) * x
  TMBad::ad_segment operator() (const std::vector<TMBad::ad_segment> &args) {
    // Unpack input
    vec A = args[0];
    vec x = args[1];
    vec N_= args[2];
    int N = (int) N_[0].Value();
    // Evaluate series
    vec term(x), y(x);
    for (int n=1; n<N; n++) {
      term = multiply(A, term) / n;
      y += term;
    }
    return y;
  }
};

template<>
struct expm_series<double> {
  typedef vectorize::vector<double> vec;
  int N;
  Eigen::SparseMatrix<double> A;
  config<double> cfg;
  void update(Eigen::SparseMatrix<double> &A) {
    this->A = A;
  }
  expm_series() {}
  expm_series(const Eigen::SparseMatrix<double> &A, double N, config<double> cfg=config<double>()) : N(N), A(A), cfg(cfg)
  {
    if ((int) N > cfg.Nmax) {
      if (cfg.warn)
        Rf_warning("expm: N terms reduced to Nmax (%i)", (int) cfg.Nmax);
      this->N = cfg.Nmax;
    }
  }
  vec operator()(vec x) {
    // Evaluate series
    vec term(x), y(x);
    for (int n=1; n<N; n++) {
      term = A * term / n;
      y += term;
    }
    return y;
  }
};

/** \brief Matrix exponential of generator matrix multiplied by vector

    This method assumes Q is a **generator matrix**, i.e
    `rowsums(Q)=0`, diagonal elements are negative, other entries are
    non-negative.  Uniformization is used to determine the number of
    terms in series.

    \note Multiplication with vector is from the left: `x^T * exp(Q)`
*/
template<class T>
struct expm_generator {
  typedef vectorize::vector<T> vec;
  expm_series<T> ExpA;
  T rho;
  /** \brief Uniformization Grassmann (1977) eq. 10. => error < 1e-4
      `rho=abs(max(diag(Q)))`
  */
  T getN(T rho) {
    return ceil(rho + 4*sqrt(rho) + 5);
  }
  // FIXME: Calculate rho=abs(max(diag(Q)))
  expm_generator(Eigen::SparseMatrix<T> &Q, config<T> cfg=config<T>()) {
    using TMBad::min;
    Eigen::SparseMatrix<T> A(Q);
    vector<T> d = Q.diagonal();
    if (d.size() > 0) {
      T M = d[0];
      for (int i=1; i<d.size(); i++) {
        M = min(M, d[i]);
      }
      rho = -M;
    }
    A.diagonal().array() += rho;
    ExpA = expm_series<T>(A, getN(rho), cfg);
  }
  /** \brief Evaluate `x^T*exp(Q)` */
  vec operator()(vec x) {
    vec y = ExpA(x);
    y = exp(-rho) * y;
    if (ExpA.cfg.normalize) {
      y = y / sum(y);
    }
    return y;
  }
};

} // End namespace

#endif
