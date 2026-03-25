#include <set>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <memory>
# ifdef _OPENMP
#include <omp.h>
# endif
#include <R.h>
#include <Rinternals.h>

/* This file contains:

   cs_ichol()              : Incomplete uplooking Cholesky and symbolic analysis.
   cs_ichol_update()       : Incomplete uplooking Cholesky factorize.
   cs_ichol_adjoint_update : Adjoint code

   cs_ichol_update_omp()         : Parallel version
   cs_ichol_adjoint_update_omp() : Parallel version (not done yet)
*/

/* --- primary CSparse routines and data structures ------------------------- */
#define csi int
typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
{
    csi nzmax ;     /* maximum number of entries */
    csi m ;         /* number of rows */
    csi n ;         /* number of columns */
    csi *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
    csi *i ;        /* row indices, size nzmax */
    double *x ;     /* numerical values, size nzmax */
    csi nz ;        /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;
/* --- 'cs_sparse' with allocated data -------------------------------------- */
/* Note: Intended to be wrapped in std::unique_ptr. Unsafe on its own! */
struct cs_sparse_alloc : cs {
  struct data {
    std::vector<csi> p;
    std::vector<csi> i;
    std::vector<double> x;
    data(size_t psize, size_t isize, size_t xsize) :
      p(psize), i(isize), x(xsize) {}
  } alloc;
  void sync() {
    this->p = alloc.p.size() ? alloc.p.data() : NULL;
    this->i = alloc.i.size() ? alloc.i.data() : NULL;
    this->x = alloc.x.size() ? alloc.x.data() : NULL;
  }
  cs_sparse_alloc(csi m, csi n, csi nzmax, csi values, csi triplet) :
    cs({std::max (nzmax, 1), m, n, NULL, NULL, NULL, triplet ? 0 : -1}),
    alloc(data(triplet ? nzmax : n+1, // p.size()
               nzmax,                 // i.size()
               values ? nzmax : 0     // x.size()
               ))
  {
    sync();
  }
};
// Replacement of 'cs*' with managed data. Note: get() member returns 'cs*'.
struct cs_ptr : std::unique_ptr<cs_sparse_alloc> {
  typedef std::unique_ptr<cs_sparse_alloc> Base;
  cs_ptr(csi m, csi n, csi nzmax, csi values, csi triplet) :
    Base(std::make_unique<cs_sparse_alloc> (m,n,nzmax,values,triplet))
  { }
};

cs_ptr cs_transpose (const cs *A, csi values) ;

/* Custom container to hold a sorted sequence of integers, optimized
   for insertions where new elements are very likely:
   - In the sequence already or
   - Greater than, and close to, previously inserted element
*/
struct sorted_ints {
  static const int NA = -1;
  int beg, prv_insert;
  // i in sequence <=> nxt[i] != NA
  std::vector<int> nxt;
  sorted_ints() {}
  sorted_ints(size_t n) : beg(n), prv_insert(n), nxt(n, NA) {}
  // Search using existing element as offset
  int find_from(int find, int from) {
    for ( ; nxt[from] <= find; from = nxt[from] );
    return from;
  }
  void insert_after(int elt, int from) {
    nxt[elt] = nxt[from];
    nxt[from] = elt;
  }
  void insert(int elt) {
    // Element already inserted
    if (nxt[elt] != NA) {
      prv_insert = elt;
      return;
    }
    // New initial element?
    if (elt < beg) {
      nxt[elt] = beg;
      beg = elt;
      return;
    }
    // Use previous insert as offset if possible
    int i = (elt > prv_insert && nxt[prv_insert] != NA ?
             prv_insert : beg);
    // Find element starting from offset
    i = find_from(elt, i);
    // Insert
    insert_after(elt, i);
    // Update
    prv_insert = elt;
  }
  const int* begin() { return &beg; }
  // erase: only first element!
  void erase(int i) {
    if ( i == beg ) {
      int tmp = nxt[beg];
      nxt[beg] = NA;
      beg = tmp;
    }
  }
  bool empty() { return (beg == (int) nxt.size()); }
};

/* Helpers for cs_ichol: 'row_pattern' and 'col_pattern' */
struct row_pattern {
  sorted_ints s; // Pattern
  std::vector<double> x; // Numerical values
  row_pattern(int n) : s(n), x(n, 0) { }
  bool empty() { return s.empty(); }
  int top() { return *s.begin(); }
  // Element access
  double& operator[](int i) { s.insert(i); return x[i]; }
  void erase(int i) { s.erase(i); x[i] = 0; }
};
struct entry {
  int index;
  double value;
};
typedef std::vector<entry> col_pattern;

/* cs_ichol
   ========

   Incomplete *up-looking* Cholesky based on Csparse routine
   'cs_chol'. That is, it calculates L in compressed column format one
   row at a time.

   - Only the upper triangle of input is accessed.
   - Output is a lower triangular matrix representing pattern and
     numerical values of LDL factorization: 'L-I+D'
   - Factor entries below tolerance 'tol' are dropped *except* entries that are already present in A.
   - Error tolerance is not to be understood as an error on the
     complete factorization: The error is considered one row at a time
     *assuming* that previous rows are without error!
   - Passing 'tol=0' gives a complete factorization.
*/
cs_ptr cs_ichol (const cs *C, double tol)
{
    double d, lki, *Lx, *Cx ;
    csi i, p, k, n, *Li, *Lp, *Cp, *Ci ;
    n = C->n ;
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    row_pattern x(n);
    std::vector<col_pattern > Lcol(n);
    std::vector<double> Ldiag(n);
    std::vector<bool> must_keep(n, false);
    for (k = 0 ; k < n ; k++)       /* compute L(k,:) for L*L' = C */
    {
        /* --- Nonzero pattern of L(k,:) ------------------------------------ */
        for (p = Cp [k] ; Ci [p] < k ; p++)       /* x = full(triu(C(:,k))) */
        {
          x [Ci [p]] = Cx [p] ;
          must_keep[Ci[p]] = true;
        }
        if (Ci[p] != k) Rf_error("ichol: rows must be sorted with diagonal present");
        d = Cx [p] ;                     /* d = C(k,k) */
        /* --- Triangular solve --------------------------------------------- */
        while (!x.empty())          /* solve L(0:k-1,0:k-1) * x = C(:,k) */
        {
            i = x.top() ;               /* s [top..n-1] is pattern of L(k,:) */
            double xi = x[i];
            x.erase (i) ;                     /* clear x for k+1st iteration */
            lki = xi / Ldiag[i];
            if (!must_keep[i] && std::abs(lki) < tol) continue;
            must_keep[i] = false; // clear
            for (auto row : Lcol[i] )
            {
                x [row.index] -= row.value * xi ;
            }
            d -= lki * lki * Ldiag[i] ;            /* d = d - L(k,i)*L(k,i) */
            Lcol[i].push_back({k, lki});   /* store L(k,i) in column i */
        }
        /* --- Compute L(k,k) ----------------------------------------------- */
        Ldiag[k] = d ;   /* store L(k,k) = sqrt (d) in column k */
    }
    /* finalize L */
    int nzmax=0;
    for (auto x : Lcol) nzmax += x.size() + 1;
    cs_ptr L (n, n, nzmax, 1, 0) ;    /* allocate result */
    if (!L) Rf_error("Failed to allocate output");
    Lp = L->p ; Li = L->i ; Lx = L->x ;
    int col=0; p=0;
    for (auto x : Lcol) {
      Lp[col] = p;
      Li[p] = col;
      Lx[p] = Ldiag[col];
      p++;
      for (auto row : x) {
        Li[p] = row.index;
        Lx[p] = row.value;
        p++;
      }
      col++;
    }
    Lp[col] = nzmax;
    /* done */
    return L;
}

/* cs_ichol_update:
   ================

   - As with 'cs_ichol', the input A is only accessed via its upper triangle!
   - The input L is lower triangular matrix obtained as output from 'cs_ichol'.
   - On output the numerical values of L are updated to hold the
     incomplete factorization of the A matrix restricted to the pattern
     of L. That is
     - L*D*L' = A on the pattern of L+L'.
     - Error can only occur on the 'residual pattern' (Pattern(L*L') \ Pattern(L+L')).
   - We define the max 'recursive error' as the max error on
     abs(L[i,]) assuming all previous rows are without error.
     On request, this error is calulated via the 'err' argument.

     Note: 'cs_ichol_update' needs the row pattern (R = Transposed of
     L) without numerical values (R -> x). This can in principle be
     precomputed, but for now we skip this optimization and calculate
     R via cs_transpose inside the function.
*/
bool cs_ichol_update (const cs *A, cs *L, double* err = NULL)
{
    double d, lki, *Lx, *Cx ;
    csi top, i, p, k, n, *Li, *Lp, *Ri, *Rp, *Cp, *Ci ;
    cs *C;
    n = A->n ;
    std::vector<csi> c(n);
    std::vector<double> x(n);
    C = (cs*) A;
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    Lp = L->p ; Li = L->i ; Lx = L->x ;
    cs_ptr R = cs_transpose(L, 0); // Pattern only
    Rp = R->p ; Ri = R->i ;
    std::vector<int> mark(n, 0);
    std::vector<int> resid;
    for (k = 0 ; k < n ; k++) c [k] = Lp [k];
    for (k = 0 ; k < n ; k++)       /* compute L(k,:) for L*L' = C */
    {
        /* --- Nonzero pattern of L(k,:) ------------------------------------ */
        x [k] = 0 ;                                 /* x (0:k) is now zero */
        for (p = Cp [k] ; p < Cp [k+1] ; p++)       /* x = full(triu(C(:,k))) */
        {
            if (Ci [p] <= k) x [Ci [p]] = Cx [p] ;
        }
        d = x [k] ;                     /* d = C(k,k) */
        x [k] = 0 ;                     /* clear x for k+1st iteration */
        /* --- Triangular solve --------------------------------------------- */
        for ( top=Rp[k]; top<Rp[k+1]; top++) mark[Ri[top]] = true;
        /* solve L(0:k-1,0:k-1) * x = C(:,k) */
        for ( top=Rp[k]; top<Rp[k+1]-1; top++)
        {
            i = Ri [top] ;               /* s [top..n-1] is pattern of L(k,:) */
            double xi = x[i];
            x [i] = 0 ;                 /* clear x for k+1st iteration */
            lki = xi / Lx [Lp [i]] ;
            if (err) {
              for (p = Lp [i] + 1 ; p < c [i] ; p++)
                {
                  x [Li [p]] -= Lx [p] * xi;
                  if (!mark[Li[p]]) {
                    resid.push_back(Li[p]);
                    mark[Li[p]] = true;
                  }
                }
            } else {
              for (p = Lp [i] + 1 ; p < c [i] ; p++) {
                x [Li [p]] -= Lx [p] * xi * mark[Li[p]];
              }
            }
            d -= lki * lki * Lx [Lp [i]] ;            /* d = d - L(k,i)*L(k,i) */
            p = c [i]++ ;
            Lx [p] = lki ;                            /* store L(k,i) in column i */
        }
        for ( top=Rp[k]; top<Rp[k+1]; top++) mark[Ri[top]] = false;
        if (err) {
          for (size_t i=0; i<resid.size(); i++) {
            *err = std::max(*err, std::abs( x[resid[i]]  / Lx[Lp[resid[i]]] ));
            x[resid[i]] = 0;
            mark[resid[i]] = false;
          }
          resid.resize(0);
        }
        /* --- Compute L(k,k) ----------------------------------------------- */
        p = c [k]++ ;
        Lx [p] = d ;                  /* store d in L(k,k) */
    }
    return true;
}

/*
  L : On input the LDL factor. On output derivative of
  f : A -> L(A) -> log(det(L(A)))

  - We loop over rows in reverse k=n-1,...,0
  - The A matrix adjoint is updated without using A itself.
  - One-by-one rows of L(A) are replaced by A adjoints.
*/
void cs_ichol_adjoint_update (cs *L) {
  double *Lx, *Cx ;
  csi top, i, p, k, n, *Li, *Lp, *Ri, *Rp, *Cp, *Ci ;
  cs *C;
  n = L->n ;
  std::vector<int> c(n);
  // Adjoints
  std::vector<double> dx(n);
  std::vector<double> dL(L->nzmax); // pattern as L
  Lp = L->p ; Li = L->i ; Lx = L->x ;
  cs_ptr R = cs_transpose(L, 0); // Pattern only
  Rp = R->p ; Ri = R->i ;
  std::vector<int> mark(n, 0);
  for (k = 0 ; k < n ; k++) c [k] = Lp [k+1];
  // Seed adjoint (derivative of log(.))
  for (k = 0 ; k < n ; k++) {
    dL[Lp[k]] += 1. / Lx[Lp[k]];
  }
  // LDL adjoint
  for (k = n ; k > 0 ; ) {
    k--; c[k]--; p = c[k];
    // 'd' and its adjoint
    double d = Lx[p]; double dd = dL[p];
    for (top = Rp[k]; top < Rp[k+1]; top++)
      mark[Ri[top]] = true;
    for (top = Rp[k+1]-1; top > Rp[k]; ) {
      top--;
      i = Ri [top] ;
      c[i]--; p = c[i];
      /* adjoint of:
         d -= lki * lki * Lx [Lp [i]]
         where d=Lx[p] has adjoint dL[p]
      */
      double lki = Lx[p]; double& dlki = dL[p];
      dlki += -2. * lki * Lx [Lp [i]] * dd; // lki adjoint
      dL[Lp[i]] -= lki * lki * dd; // Lx [Lp [i]] adjoint
      double xi = lki * Lx[Lp[i]]; // restore
      double dxi = 0;
      for (p = c [i]; p > Lp [i] + 1 ; ) {
        p--;
        /* x [Li [p]] -= Lx [p] * xi * mark[Li[p]]; */
        dxi -= Lx[p] * dx[Li[p]] * mark[Li[p]];
        dL[p] -= dx[Li[p]] * xi * mark[Li[p]];
      }
      // lki = xi / Lx [Lp [i]] ;
      dxi += dlki / Lx [Lp [i]] ;
      dL[Lp[i]] -= dlki * xi / (Lx [Lp [i]]*Lx [Lp [i]]);
      // double xi = x[i];
      dx[i] += dxi;
    }
    for ( top=Rp[k]; top<Rp[k+1]; top++) mark[Ri[top]] = false;
    // d = x[k]
    dx[k] += dd;
    for (p = Rp [k] ; p < Rp [k+1] ; p++) {
      // dA adjoint: store in k'th row af dL which is not needed anymore
      int i = Ri[p];
      if (i <= k) {
        dL[c[i]] = dx[i] ;
      }
      dx[i] = 0; // clear
    }
  }
  Memcpy(L->x, (double*)(dL.data()), dL.size());
}

/* --- Parallel versions ---------------------------------------------------- */

struct thread_data {
  std::vector<int> mark;
  std::vector<int> resid;
  std::vector<double> x;
  thread_data() {}
  thread_data(size_t n) : mark(n, 0), x(n, 0) {}
};
void cs_ichol_update_omp_wrk (const cs *A, cs *L, double* err,
                              const cs *R, int K, int N,
                              std::vector<int> &c,
                              thread_data &tdat) {
  std::vector<double> &x = tdat.x;
  std::vector<int> &mark = tdat.mark;
  std::vector<int> &resid = tdat.resid;
  double d, lki, *Lx, *Cx ;
  csi top, i, p, k, n, *Li, *Lp, *Ri, *Rp, *Cp, *Ci ;
  cs *C;
  C = (cs*) A;
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  Lp = L->p ; Li = L->i ; Lx = L->x ;
  Rp = R->p ; Ri = R->i ;
  for (k = K ; k < K+N ; k++)       /* compute L(k,:) for L*L' = C */
    {
      /* --- Nonzero pattern of L(k,:) ------------------------------------ */
        x [k] = 0 ;                                 /* x (0:k) is now zero */
        for (p = Cp [k] ; p < Cp [k+1] ; p++)       /* x = full(triu(C(:,k))) */
        {
            if (Ci [p] <= k) x [Ci [p]] = Cx [p] ;
        }
        d = x [k] ;                     /* d = C(k,k) */
        x [k] = 0 ;                     /* clear x for k+1st iteration */
        /* --- Triangular solve --------------------------------------------- */
        for ( top=Rp[k]; top<Rp[k+1]; top++) mark[Ri[top]] = true;
        /* solve L(0:k-1,0:k-1) * x = C(:,k) */
        for ( top=Rp[k]; top<Rp[k+1]-1; top++)
        {
            i = Ri [top] ;               /* s [top..n-1] is pattern of L(k,:) */
            double xi = x[i];
            x [i] = 0 ;                 /* clear x for k+1st iteration */
            lki = xi / Lx [Lp [i]] ;
            if (err) {
              for (p = Lp [i] + 1 ; p < c [i] ; p++)
                {
                  x [Li [p]] -= Lx [p] * xi;
                  if (!mark[Li[p]]) {
                    resid.push_back(Li[p]);
                    mark[Li[p]] = true;
                  }
                }
            } else {
              for (p = Lp [i] + 1 ; p < c [i] ; p++) {
                x [Li [p]] -= Lx [p] * xi * mark[Li[p]];
              }
            }
            d -= lki * lki * Lx [Lp [i]] ;            /* d = d - L(k,i)*L(k,i) */
            p = c [i]++ ;
            Lx [p] = lki ;
        }
        for ( top=Rp[k]; top<Rp[k+1]; top++) mark[Ri[top]] = false;
        if (err) {
          for (size_t i=0; i<resid.size(); i++) {
            *err = std::max(*err, std::abs( x[resid[i]]  / Lx[Lp[resid[i]]] ));
            x[resid[i]] = 0;
            mark[resid[i]] = false;
          }
          resid.resize(0);
        }
        /* --- Compute L(k,k) ----------------------------------------------- */
        p = c [k]++ ;
        Lx [p] = d ;
    }
}

cs_ptr cs_parallel_schedule (const cs *L) { // Input R=t(L)
  double d, lki, *Lx, *Cx ;
  csi top, i, p, k, n, *Li, *Lp, *Ri, *Rp, *Cp, *Ci ;
  cs *C;
  n = L->n ;
  Lp = L->p ; Li = L->i ; Lx = L->x ;
  cs_ptr R = cs_transpose(L, 0); // Pattern only
  Rp = R->p ; Ri = R->i ;
  static const int NA = -1;
  std::vector<int> mark(n, NA);
  std::vector<int> stack; // = which(mark != -1)
  int cur_max = NA;
  // result
  std::vector<int> rows;
  std::vector<int> chunk(1, 0);
  for (int k = 0; k<n; k++) {
    bool all_NA = true;
    bool all_max = true;
    for (p = Rp[k]; p < Rp[k+1]; p++) {
      int i = Ri[p];
      if (mark[i] != NA) {
        all_NA = false;
        if (mark[i] != cur_max)
          all_max = false;
      }
    }
    if (all_NA) {
      // Start new row
      cur_max = k;
      rows.push_back(k);
    } else if (!all_max) {
      // conflict
      // Start new row *and* new chunk
      chunk.push_back(rows.size()); // New chunk start here
      cur_max = k;
      rows.push_back(k);
      for (size_t i = 0; i < stack.size(); i++) mark[stack[i]] = NA;
      stack.resize(0);
    }
    // Mark kth row
    for (p = Rp[k]; p < Rp[k+1]; p++) {
      int i = Ri[p];
      if (mark[i] == NA) stack.push_back(i);
      mark[i] = cur_max;
    }
  }
  chunk.push_back(rows.size()); // Must end like this
  cs_ptr S(n, chunk.size() - 1, rows.size(), 0, 0) ;
  Memcpy(S->i, (int*)(rows.data()), rows.size());
  Memcpy(S->p, (int*)(chunk.data()), chunk.size());
  return S;
}

bool cs_ichol_update_omp (const cs *A, cs *L, cs *S, double* err = NULL)
{
    double d, lki, *Lx, *Cx ;
    csi top, i, p, k, n, *Li, *Lp, *Ri, *Rp, *Cp, *Ci, *Sp, *Si ;
    cs *C;
    n = A->n ;
    std::vector<csi> c(n);
    C = (cs*) A;
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    Lp = L->p ; Li = L->i ; Lx = L->x ;
    cs_ptr R = cs_transpose(L, 0); // Pattern only
    Rp = R->p ; Ri = R->i ;
    Sp = S->p ; Si = S->i ;
    for (k = 0 ; k < n ; k++) c [k] = Lp [k];
    // Thread private data
    std::vector<thread_data> tdat(omp_get_max_threads());
#pragma omp parallel
    {
      // Initialize thread data in parallel (first touch)
      tdat[omp_get_thread_num()] = thread_data(n);
    }
    for (int batch = 0; batch < S->n; batch++) {
#pragma omp parallel for
      for (p = Sp[batch] ; p < Sp[batch+1]; p++) {
        int k1 = Si[p];
        bool last_p = (p == Sp[S->n] - 1);
        int k2 = ( last_p ? n : Si[p+1]);
        cs_ichol_update_omp_wrk (A, L, err,
                                 R.get(),
                                 k1, k2-k1,
                                 c,
                                 tdat[omp_get_thread_num()]);
      }
    }
    return true;
}



/* Right-looking incomplete Cholesky factorization
   ===============================================

   - Input: Lower triangular matrix representing tril(Q) and lower
     triangular pattern matrix L (output of 'cs_ichol').
   - Output: Numerical values of L are updated with LDL factorization.
   - Note: (L,D) repesented in a single matrix as 'L-I+D'.
*/
void cs_chol_right(cs* Q, cs* L) {
  // Pointers
  int* Qp = (int*) Q->p;
  int* Qi = (int*) Q->i;
  double* Qx = (double*) Q->x;
  int* Lp = (int*) L->p;
  int* Li = (int*) L->i;
  double* Lx = (double*) L->x;
  int nc = Q->n;
  std::vector<double> buf(nc);
  double *w = buf.data();
  // Embed Q in L
  for (int col = 0; col < nc; col++) {
    // Scatter
    for (int p = Qp[col]; p < Qp[col+1]; p++) w[Qi[p]] = Qx[p];
    // Gather
    for (int p = Lp[col]; p < Lp[col+1]; p++) Lx[p] = w[Li[p]];
    // Clear
    for (int p = Qp[col]; p < Qp[col+1]; p++) w[Qi[p]] = 0;
  }
  // Right looking Cholesky
  for (int col = 0; col < nc; col++) {
    int p_diag = Lp[col];
    double d = Lx[p_diag];
    // Compute and scatter this off-diagonal column
    for (int p = Lp[col] + 1; p < Lp[col+1]; p++) {
      Lx[p] /= d;
      w[Li[p]] = Lx[p];
    }
    // Loop through target subset (non-zeros of relevant columns)
    for (int p = Lp[col] + 1; p < Lp[col+1]; p++) {
      int tcol = Li[p]; // target column
      for (int tp = Lp[tcol]; tp < Lp[tcol+1]; tp++) {
        int trow = Li[tp]; // target row
        Lx[tp] -= w[trow] * w[tcol] * d ;
      }
    }
    // Reset workspace
    for (int p = Lp[col] + 1; p < Lp[col+1]; p++) {
      w[Li[p]] = 0;
    }
  }
}

/* Transform LDL to the derivative: d/dQ(sum(log(diag(D))))

   This function calculates the adjoint code of 'cs_chol_right(Q, L)'
   in the range direction given by D^-1, i.e. the derivative of
   d/dD(sum(log(diag(D)))).

   Precomputed LDL factor (L-I+D) is supplied as input, and
   overwritten with derivatives on output.
*/
void cs_dchol_right(cs* L) {
  // Pointers
  int* Lp = (int*) L->p;
  int* Li = (int*) L->i;
  double* Lx = (double*) L->x;
  int nc = L->n;
  std::vector<double> buf(nc);
  double *w = buf.data();
  std::vector<double> dbuf(nc);
  double *dw = dbuf.data();
  std::vector<int> mark(nc, 0);
  for (int col = nc; col > 0; ) {
    col--;
    // Get d
    int p_diag = Lp[col];
    double d = Lx[p_diag];
    // scatter this off-diagonal column
    for (int p = Lp[col] + 1; p < Lp[col+1]; p++) {
      w[Li[p]] = Lx[p];
      mark[Li[p]] = true;
    }
    // Loop through target subset (non-zeros of relevant columns)
    for (int p = Lp[col] + 1; p < Lp[col+1]; p++) {
      int tcol = Li[p]; // target column
      for (int tp = Lp[tcol]; tp < Lp[tcol+1]; tp++) {
        int trow = Li[tp]; // target row
        if (mark[trow] && mark[tcol]) {
          dw[trow] -= Lx[tp] * w[tcol];
          dw[tcol] -= Lx[tp] * w[trow];
          dw[col] += Lx[tp] * w[trow] * w[tcol];
        }
      }
    }
    dw[col] += 1. / d;
    // Gather this column
    for (int p = Lp[col] ; p < Lp[col+1]; p++) {
      Lx[p] = dw[Li[p]];
    }
    // Reset workspaces
    for (int p = Lp[col] ; p < Lp[col+1]; p++) {
      w[Li[p]] = 0;
      dw[Li[p]] = 0;
      mark[Li[p]] = false;
    }
  }
}

/* solve Lx=b where x and b are dense.  x=b on input, solution on output. */
template<bool unit=false, bool diag=false>
csi cs_lsolve (const cs *L, double *x)
{
    csi p, j, n, *Lp, *Li ;
    double *Lx ;
    n = L->n ; Lp = L->p ; Li = L->i ; Lx = L->x ;
    for (j = 0 ; j < n ; j++)
    {
        if (!unit) x [j] /= Lx [Lp [j]] ;
        if (!diag) {
          for (p = Lp [j]+1 ; p < Lp [j+1] ; p++)
            {
              x [Li [p]] -= Lx [p] * x [j] ;
            }
        }
    }
    return (1) ;
}
/* solve L'x=b where x and b are dense.  x=b on input, solution on output. */
template<bool unit=false, bool diag=false>
csi cs_ltsolve (const cs *L, double *x)
{
    csi p, j, n, *Lp, *Li ;
    double *Lx ;
    n = L->n ; Lp = L->p ; Li = L->i ; Lx = L->x ;
    for (j = n-1 ; j >= 0 ; j--)
    {
        if (!diag) {
          for (p = Lp [j]+1 ; p < Lp [j+1] ; p++)
            {
              x [j] -= Lx [p] * x [Li [p]] ;
            }
        }
        if (!unit) x [j] /= Lx [Lp [j]] ;
    }
    return (1) ;
}
/* solve LDL'x=b */
void cs_solve (const cs *L, double *x) {
  cs_lsolve  <1,0> (L, x);  // L^-1 x
  cs_lsolve  <0,1> (L, x);  // D^-1 x
  cs_ltsolve <1,0> (L, x);  // L'^-1 x
}

/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
double cs_cumsum (csi *p, csi *c, csi n)
{
    csi i, nz = 0 ;
    double nz2 = 0 ;
    if (!p || !c) return (-1) ;     /* check inputs */
    for (i = 0 ; i < n ; i++)
    {
        p [i] = nz ;
        nz += c [i] ;
        nz2 += c [i] ;              /* also in double to avoid csi overflow */
        c [i] = p [i] ;             /* also copy p[0..n-1] back into c[0..n-1]*/
    }
    p [n] = nz ;
    return (nz2) ;                  /* return sum (c [0..n-1]) */
}

/* C = A' */
cs_ptr cs_transpose (const cs *A, csi values)
{
    csi p, q, j, *Cp, *Ci, n, m, *Ap, *Ai ;
    double *Cx, *Ax ;
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    cs_ptr C(n, m, Ap [n], values && Ax, 0) ;              /* allocate result */
    std::vector<csi> w(m);                                 /* get workspace */
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (p = 0 ; p < Ap [n] ; p++) w [Ai [p]]++ ;          /* row counts */
    cs_cumsum (Cp, w.data(), m) ;                          /* row pointers */
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            Ci [q = w [Ai [p]]++] = j ; /* place A(i,j) as entry C(j,i) */
            if (Cx) Cx [q] = Ax [p] ;
        }
    }
    return C;
}

// R interface
static R_INLINE
SEXP ALLOC_SLOT(SEXP obj, SEXP nm, SEXPTYPE type, R_xlen_t length)
{
  SEXP val = allocVector(type, length);
  R_do_slot_assign(obj, nm, val);
  return val;
}
void cs2r(SEXP ans, cs* L) {
  /* Dim : */
  int* Dim = INTEGER(ALLOC_SLOT(ans, Rf_install("Dim"), INTSXP, 2));
  Dim[0] = L->m; Dim[1] = L->n;
  /* p : */
  int* p = INTEGER(ALLOC_SLOT(ans, Rf_install("p"), INTSXP, L->n + 1));
  Memcpy(p, (int*)(L->p), L->n + 1);
  /* i : */
  int* i = INTEGER(ALLOC_SLOT(ans, Rf_install("i"), INTSXP, L->nzmax));
  Memcpy(i, (int*)(L->i), L->nzmax);
  /* x : */
  if (L->x != NULL) {
    double* x = REAL(ALLOC_SLOT(ans, Rf_install("x"), REALSXP, L->nzmax));
    Memcpy(x, (int*)(L->x), L->nzmax);
  }
}
cs r2cs(SEXP X) {
  cs A;
  A.nzmax = LENGTH(R_do_slot(X, Rf_install("i")));
  A.m = INTEGER(R_do_slot(X, Rf_install("Dim")))[0];
  A.n = INTEGER(R_do_slot(X, Rf_install("Dim")))[1];
  A.p = INTEGER(R_do_slot(X, Rf_install("p")));
  A.i = INTEGER(R_do_slot(X, Rf_install("i")));
  bool has_x = R_has_slot(X, Rf_install("x"));
  A.x = has_x ? REAL(R_do_slot(X, Rf_install("x"))) : NULL;
  A.nz = -1;
  return A;
}

/* --- R-interface ---------------------------------------------------------- */

extern "C"
SEXP tmb_isolve(SEXP X, SEXP Y) {
  cs A = r2cs(X);
  if (A.m != A.n || A.m != LENGTH(Y))
    Rf_error("Non-conformable arguments");
  SEXP ans = PROTECT(Rf_duplicate(Y));
  double* x = REAL(ans);
  cs_solve(&A, x);
  UNPROTECT(1);
  return ans;
}
extern "C"
SEXP tmb_ichol(SEXP X, SEXP tol) {
  cs A = r2cs(X);
  double chol_tol = REAL(tol)[0];
  cs_ptr L = cs_ichol(&A, chol_tol);
  SEXP claes = PROTECT(R_do_MAKE_CLASS("dgCMatrix"));  // OR dtCMatrix ?
  SEXP ans = PROTECT(R_do_new_object(claes));
  cs2r(ans, L.get());
  UNPROTECT(2);
  return ans;
}
// Parallel schedule
extern "C"
SEXP tmb_parallel_schedule(SEXP L) {
  cs A = r2cs(L);
  cs_ptr S = cs_parallel_schedule(&A);
  SEXP claes = PROTECT(R_do_MAKE_CLASS("ngCMatrix"));
  SEXP ans = PROTECT(R_do_new_object(claes));
  cs2r(ans, S.get());
  UNPROTECT(2);
  return ans;
}
// Update incomplete Cholesky factor L
// Attribute 'parallel_schedule' selects parallel version.
extern "C"
SEXP tmb_ichol_update(SEXP X, SEXP Y, SEXP get_error) {
  bool get_err = INTEGER(get_error)[0];
  cs A = r2cs(X);
  cs L = r2cs(Y);
  double err = 0;
  SEXP SCH = Rf_getAttrib(Y, Rf_install("parallel_schedule"));
  if (!Rf_isNull(SCH)) {
    cs S = r2cs(SCH);
    cs_ichol_update_omp(&A, &L, &S, get_err ? &err : NULL);
  } else {
    cs_ichol_update(&A, &L, get_err ? &err : NULL);
  }
  SEXP error = PROTECT(Rf_ScalarReal(err));
  Rf_setAttrib(Y, Rf_install("error"), error);
  UNPROTECT(1);
  return Y;
}
// derivatives
extern "C"
SEXP tmb_ichol_adjoint_update(SEXP L) {
  SEXP X = PROTECT(Rf_duplicate(L));
  cs A = r2cs(X);
  cs_ichol_adjoint_update(&A);
  UNPROTECT(1);
  return X;
}
// Update incomplete LDL Cholesky factor L
extern "C"
SEXP tmb_ldl_update(SEXP X, SEXP Y) {
  cs A = r2cs(X);
  cs L = r2cs(Y);
  cs_chol_right(&A, &L);
  return Y;
}
extern "C"
SEXP tmb_ldl_deriv(SEXP L) {
  SEXP X = PROTECT(Rf_duplicate(L));
  cs A = r2cs(X);
  cs_dchol_right(&A);
  UNPROTECT(1);
  return X;
}
