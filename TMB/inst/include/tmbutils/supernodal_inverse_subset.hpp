#include <Eigen/CholmodSupport>
#include <unsupported/Eigen/AutoDiff>

namespace Eigen {

/** \brief Supernodal Cholesky factor with access to protected members */
struct Accessible_CholmodSupernodalLLT :
    Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double> > {
  typedef Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double> > Base;
  // Inherit all CTORs
  using Base::Base;
  cholmod_factor* get_factor() {
    return m_cholmodFactor;
  }
  StorageIndex* get_super() {
    return static_cast<StorageIndex*>(get_factor()->super);
  }
  StorageIndex* get_pi() {
    return static_cast<StorageIndex*>(get_factor()->pi);
  }
  StorageIndex* get_s() {
    return static_cast<StorageIndex*>(get_factor()->s);
  }
  StorageIndex* get_Perm() {
    return static_cast<StorageIndex*>(get_factor()->Perm);
  }
  StorageIndex* get_px() {
    return static_cast<StorageIndex*>(get_factor()->px);
  }
  size_t get_n() {
    return get_factor()->n;
  }
  size_t get_nsuper() {
    return get_factor()->nsuper;
  }
  size_t get_xsize() {
    return get_factor()->xsize;
  }
  double* get_x() {
    return static_cast<double*>(get_factor()->x);
  }
};
template<class T>
struct SupernodalInverseSubset {
  // typedefs
  typedef Accessible_CholmodSupernodalLLT Base;
  typedef Base::StorageIndex StorageIndex;
  typedef SparseMatrix<T> SpMat;
  typedef Matrix<T, Dynamic, Dynamic> Mat;
  typedef Map<Array<StorageIndex, Dynamic, 1> > Rows_ref;
  enum Operation {
    Get,
    Set,
    Update
  };
  // workspaces
  std::shared_ptr<Base> chm;
  std::vector<StorageIndex> col2super;
  std::vector<StorageIndex> col2offset;
  std::vector<StorageIndex> idg;
  std::vector<StorageIndex> iwrk;
  std::vector<T> ans;
  // Helper functions
  void init_col2super() {
    if (col2super.size() > 0) return;
    col2super.resize( chm->get_n() );
    StorageIndex *super = chm->get_super();
    Index nsuper = chm->get_nsuper();
    for (Index k=0, l=0; k < nsuper; k++) {
      StorageIndex ncols = super[k + 1] - super[k];
      while (ncols--) col2super[l++] = k;
    }
  }
  void init_col2offset() {
    if (col2offset.size() > 0) return;
    col2offset.resize( chm->get_n(), 0 );
    StorageIndex *pi = chm->get_pi();
    for (size_t i = 1; i < col2offset.size(); i++) {
      StorageIndex k = col2super[i - 1];
      StorageIndex nrows = pi[k + 1] - pi[k];
      col2offset[i] = col2offset[i - 1] + nrows; 
    }
  }
  void chm_index_scatter(StorageIndex col, StorageIndex* wrk) {
    init_col2super();
    init_col2offset();
    StorageIndex k = col2super[col];
    StorageIndex *pi = chm->get_pi();
    StorageIndex *s = chm->get_s();
    for (StorageIndex i = pi[k], j = col2offset[col] ; i < pi[k + 1] ; i++) {
      wrk[s[i]] = j++;
    }
  }
  std::vector<StorageIndex> index_gather(const SpMat &mat) {
    std::vector<StorageIndex> ans;
    StorageIndex *Perm = chm->get_Perm();
    std::vector<StorageIndex> iPerm(mat.rows());
    for (size_t i=0; i<iPerm.size(); i++)
      iPerm[Perm[i]] = i;
    std::vector<StorageIndex> wrk(mat.rows());
    for (StorageIndex k=0; k<mat.outerSize(); ++k) {
      chm_index_scatter(iPerm[k], wrk.data());
      for (typename SpMat::InnerIterator it(mat, k); it; ++it) {
        if (iPerm[it.row()] < iPerm[it.col()]) ans.push_back(-1); // Skip upper triangle
        else ans.push_back( wrk[iPerm[it.row()]] );
      }
    }
    return ans;
  }
  template<Operation Op>
  void values(SpMat &mat) {
    if (idg.size() == 0)
      idg = index_gather(mat);
    T* vptr = mat.valuePtr();
    for (size_t i=0; i<idg.size(); i++) {
      if (idg[i] != -1) {
        if (Op == Set)
          ans[idg[i]] = vptr[i];
        if (Op == Get)
          vptr[i] = ans[idg[i]];
      }
    }
  }
  Rows_ref rws(StorageIndex k) {
    StorageIndex *pi = chm->get_pi();
    StorageIndex *s = chm->get_s();
    StorageIndex nrows = pi[k + 1] - pi[k];
    return Rows_ref(s + pi[k], nrows);
  }
  StorageIndex nrow(StorageIndex k) {
    StorageIndex *pi = chm->get_pi();
    StorageIndex nrows = pi[k + 1] - pi[k];
    return nrows;
  }
  StorageIndex ncol(StorageIndex k) {
    StorageIndex *super = chm->get_super();
    StorageIndex ncols = super[k + 1] - super[k];
    return ncols;
  }
  template<Operation Op>
  void denseBlock(Mat &ans,
                  T* x,
                  StorageIndex* r,
                  StorageIndex* c,
                  StorageIndex nr,
                  StorageIndex nc) {
    init_col2super();
    if (iwrk.size() != chm->get_n())
      iwrk.resize(chm->get_n());
    // Super node specific information
    std::vector<StorageIndex> rel_rows; // rows relative to supernode
    StorageIndex i_start = 0;
    for (StorageIndex j = 0; j < nc; j++) {
      StorageIndex k = col2super[c[j]];
      auto X = getSuperNode(x, k);
      StorageIndex *super = chm->get_super();
      StorageIndex rel_col = c[j] - super[k]; // column relative to this supernode
      if (j == 0 || k != col2super[c[j-1]]) {
        // New supernode: Must update relative rows
        // Scatter relative supernode rows
        Rows_ref rows = rws(k);
        for (StorageIndex i=0; i<rows.size(); i++)
          iwrk[rows[i]] = i;
        // Gather relative supernode rows
        rel_rows.resize(nr);
        for (StorageIndex i=0; i<nr; i++)
          rel_rows[i] = iwrk[r[i]];
        // Find i_start
        i_start = 0;
        while (i_start < nr && r[i_start] < c[j]) i_start++;
      }
      for (StorageIndex i = i_start; i < nr; i++) {
        if (Op == Get)
          ans(i, j) = X(rel_rows[i], rel_col);
        if (Op == Set)
          X(rel_rows[i], rel_col) = ans(i, j);
        if (Op == Update)
          X(rel_rows[i], rel_col) += ans(i, j);
      }
    }
  }
  Map<Mat> getSuperNode(T* x,
                        StorageIndex k) {
    StorageIndex *px = chm->get_px();
    T* data = x + px[k];
    Map<Mat> ans(data, nrow(k), ncol(k));
    return ans;
  }
  /* --- chol: get result from cholmod -------------------------------------- */
  void chol_get_cached_values() {
    double* x = chm->get_x();
    size_t n = chm->get_xsize();
    ans.resize(n);
    for (size_t i=0; i<n; i++) ans[i] = x[i];
  }
  /* --- chol: forward recursions ------------------------------------------- */
  void chol(SpMat x) {
    // Zero initialized workspace to mimic get_factor()->x
    ans.resize(0); ans.resize(chm->get_xsize(), 0);
    // Fill SpMat into workspace
    values<Set>(x);
    // Forward loop
    StorageIndex *super = chm->get_super();
    Index nsuper = chm->get_nsuper();
    for (Index k=0; k < nsuper; k++) {
      Rows_ref rows = rws(k);
      // Notation
      StorageIndex ns = super[k + 1] - super[k];
      StorageIndex np = rows.size() - ns;
      StorageIndex* p = rows.data() + ns;
      // r:=c(s,p)
      auto Xrs = getSuperNode(ans.data(), k); // Write access !
      Xrs.block(0, 0, ns, ns) =
        Xrs.block(0, 0, ns, ns).llt().matrixL();
      if (np > 0) {
        Mat Xpp(np, np);
        Xpp.setZero();
        auto L = Xrs.block(0, 0, ns, ns);
        auto Xps = Xrs.block(ns, 0, np, ns);
        Xrs.block(ns, 0, np, ns) =
          L.template triangularView<Lower>().
          solve( Xps.transpose() ).
          transpose();
        Xpp.template selfadjointView<Lower>().
          rankUpdate(Xrs.block(ns, 0, np, ns), -1.);
        denseBlock<Update>(Xpp, ans.data(), p, p, np, np);
      }
    }
  }
  /* --- chol2inv: reverse recursions --------------------------------------- */
  void chol2inv() {
    // Reverse loop
    StorageIndex *super = chm->get_super();
    Index nsuper = chm->get_nsuper();
    for (Index k=nsuper; k > 0; ) {
      k--;
      Rows_ref rows = rws(k);
      // Notation
      StorageIndex ns = super[k + 1] - super[k];
      StorageIndex np = rows.size() - ns;
      StorageIndex* p = rows.data() + ns;
      // r:=c(s,p)
      auto Xrs = getSuperNode(ans.data(), k); // Write access !
      auto Lss = Xrs.block(0, 0, ns, ns);
      // Xss = chol2inv(Lss)
      Mat Lss_inv(ns, ns);
      Lss_inv.setIdentity();
      Lss.template triangularView<Lower>().
        solveInPlace(Lss_inv);
      Xrs.block(0, 0, ns, ns).setZero();
      Xrs.block(0, 0, ns, ns).
        template selfadjointView<Lower>().
        rankUpdate(Lss_inv.transpose(), 1.);
      if (np > 0) {
        Mat Xpp(np, np);
        denseBlock<Get>(Xpp, ans.data(), p, p, np, np);
        Mat Lps = Xrs.block(ns, 0, np, ns);
        Mat tmp =
          Lss_inv.template triangularView<Lower>().
          transpose() * Lps.transpose();
        Xrs.block(ns, 0, np, ns) =
          (-tmp * Xpp.template selfadjointView<Lower>()).transpose();
        Xrs.block(0, 0, ns, ns) -= tmp * Xrs.block(ns, 0, np, ns);
      }
    }
  }
  SpMat operator()(SpMat x) {
    // forward recursions (use cached values if T=double)
    if (isDouble<T>::value) {
      chol_get_cached_values();
    } else {
      chol(x);
    }
    // Reverse
    chol2inv();
    // Get result
    x = x * 0;
    values<Get> (x);
    return x;
  }
  SupernodalInverseSubset(std::shared_ptr<Base> chm) :
    chm(chm) { }
  void print_common() {
    cholmod_print_common("C", &(chm->cholmod()));
  }
};
} // End namespace Eigen
