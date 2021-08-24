namespace Eigen {

template<class T>
struct SimplicialInverseSubset {
  // typedefs
  typedef SimplicialLLT<SparseMatrix<T> > Base;
  typedef SimplicialLLT<SparseMatrix<double> > BaseD;
  typedef typename Base::StorageIndex StorageIndex;
  typedef typename Base::Scalar Scalar;
  typedef SparseMatrix<T> SpMat;
  enum Operation {
    //IndexScatter,
    Scatter,
    Zero,
    InnerProduct
  };
  // workspaces
  std::shared_ptr<Base> factor;
  std::vector<StorageIndex> idg;
  std::vector<StorageIndex> index_gather(const SpMat &mat) {
    SpMat L = factor->matrixL();
    const StorageIndex* Lp = L.outerIndexPtr();
    const StorageIndex* Li = L.innerIndexPtr();
    std::vector<StorageIndex> ans;
    const StorageIndex *Perm = factor->permutationP().indices().data();
    std::vector<StorageIndex> wrk(mat.rows());
    for (StorageIndex k=0; k<mat.outerSize(); ++k) {
      // Index scatter
      StorageIndex j = Perm[k];
      for (StorageIndex l=Lp[j]; l<Lp[j+1]; l++) {
        StorageIndex i = Li[l];
        wrk[i] = l;
      }
      for (typename SpMat::InnerIterator it(mat, k); it; ++it) {
        if (Perm[it.row()] < Perm[it.col()]) ans.push_back(-1); // Skip upper triangle
        else ans.push_back( wrk[Perm[it.row()]] );
      }
    }
    return ans;
  }
  void valuesGet(SpMat &mat, const SpMat &S) {
    if (idg.size() == 0)
      idg = index_gather(mat);
    T* vptr = mat.valuePtr();
    const T* ans = S.valuePtr();
    for (size_t i=0; i<idg.size(); i++) {
      if (idg[i] != -1) {
        vptr[i] = ans[idg[i]];
      }
    }
  }
  template<Operation Op, class Scalar_, class Scalar>
  Scalar column(SparseMatrix<Scalar_> &L,
                StorageIndex j,
                Scalar* x) {
    const StorageIndex* Lp = L.outerIndexPtr();
    const StorageIndex* Li = L.innerIndexPtr();
    Scalar_* Lx = L.valuePtr();
    Scalar s = 0;
    for (StorageIndex k=Lp[j]; k<Lp[j+1]; k++) {
      StorageIndex i = Li[k];
      // if (Op == IndexScatter) {
      //   x[i] = k;
      // }
      if (Op == Scatter) {
        x[i] = Lx[k];
      }
      if (Op == Zero) {
        x[i] = 0.;
      }
      if (Op == InnerProduct) {
        s += Lx[k] * x[i];
      }
    }
    return s;
  }
  Eigen::SparseMatrix<StorageIndex> LT;
  void init_transpose(SpMat L) {
    if (LT.rows() > 0) return;
    std::vector<StorageIndex> x(L.nonZeros());
    for (size_t i=0; i<x.size(); i++) x[i] = i;
    Eigen::SparseMatrix<StorageIndex> TMP =
      Eigen::Map<const Eigen::SparseMatrix<StorageIndex> > (L.rows(),
                                                            L.cols(),
                                                            L.nonZeros(),
                                                            L.outerIndexPtr(),
                                                            L.innerIndexPtr(),
                                                            x.data(),
                                                            L.innerNonZeroPtr());
    LT = TMP.transpose();
  }
  SpMat chol2inv() {
    SpMat L = factor->matrixL();
    init_transpose(L);
    // Allocate result
    SpMat S(L);
    for (StorageIndex i=0; i<S.nonZeros(); i++)
      S.valuePtr()[i] = 0;
    const StorageIndex ncol = L.cols();
    // L pointers
    const StorageIndex* Lp = L.outerIndexPtr();
    Scalar* Lx = L.valuePtr();
    // S pointers (pattern same as L)
    Scalar* Sx = S.valuePtr();
    // LT pointers
    const StorageIndex* LTp = LT.outerIndexPtr();
    StorageIndex* LTi = LT.innerIndexPtr();
    StorageIndex* LTx = LT.valuePtr();
    // Workspace dense row
    std::vector<T> wrk(ncol, 0);
    T* S_row = wrk.data();
    // Recursions
    for (StorageIndex j = ncol; j > 0; ) {
      j--;
      // Scatter sub-column S(j:n,j)
      column<Scatter> (S, j, S_row);
      // Diagonal element S(j,j)
      Scalar Sjj = 0;
      for (StorageIndex k = Lp[j] + 1; k < Lp[j+1]; k++) {
        Sjj += Lx[k] * Sx[k];
      }
      Scalar Ljj_inv = 1. / Lx[Lp[j]];
      S_row[j] = -Ljj_inv * Sjj + Ljj_inv * Ljj_inv;
      // Recursions for j'th row S(j,1:(j-1))
      for (StorageIndex k = LTp[j+1]-1; k > LTp[j]; ) {
        k--;
        StorageIndex i = LTi[k];
        S_row[i] = (-1. / Lx[Lp[i]]) * column<InnerProduct> (L, i, S_row);
      }
      for (StorageIndex k = LTp[j]; k < LTp[j+1]; k++) {
        // Insert result in sparse output
        StorageIndex i = LTi[k];
        StorageIndex kT = LTx[k];
        Sx[kT] = S_row[i];
      }
      // Clear workspace
      column<Zero> (L,  j, S_row);
      column<Zero> (LT, j, S_row);
    }
    return S;
  }
  SpMat operator()(SpMat x) {
    // Factor initialize
    if (!factor) {
      factor = std::make_shared<Base>(x);
    }
    // forward recursions
    factor->factorize(x);
    // Reverse
    SpMat S = chol2inv();
    // Get result
    x = x * 0;
    valuesGet (x, S);
    return x;
  }
  // No simple way to change factor<double> to other numerical type
  // SimplicialInverseSubset(std::shared_ptr<Base> factor) :
  //   factor(factor) { }
  SimplicialInverseSubset(std::shared_ptr<BaseD> factor) { }
};

} // Eigen
