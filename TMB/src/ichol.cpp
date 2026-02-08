#include <iostream>
#include <set>
#include <vector>
#include <cstdlib>
#include <R.h>
#include <Rinternals.h>
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
#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
cs *cs_spalloc (csi m, csi n, csi nzmax, csi values, csi triplet) ;
cs *cs_spfree (cs *A) ;
cs *cs_transpose (const cs *A, csi values) ;
cs *cs_done (cs *C, void *w, void *x, csi ok) ;
/* wrapper for malloc */
void *cs_malloc (csi n, size_t size)
{
    return (malloc (CS_MAX (n,1) * size)) ;
}
/* wrapper for calloc */
void *cs_calloc (csi n, size_t size)
{
    return (calloc (CS_MAX (n,1), size)) ;
}
/* wrapper for free */
void *cs_free (void *p)
{
    if (p) free (p) ;       /* free p if it is not already NULL */
    return (NULL) ;         /* return NULL to simplify the use of cs_free */
}

#define SIZEOF(type) type
#define CS_CALLOC(n, type) (type*) cs_calloc (n, sizeof (type))
#define CS_MALLOC(n, type) (type*) cs_malloc (n, sizeof (type))
#define CS_SPALLOC(...) (cs*) cs_spalloc ( __VA_ARGS__ )
#define CS_SPREALLOC(...) cs_sprealloc ( __VA_ARGS__ )
#define CS_DALLOC(...) cs_dalloc ( __VA_ARGS__ )
#define CS_REALLOC(X1,X2,type,X3) (type*) cs_realloc (X1, X2, sizeof(type), X3)

struct row_pattern {
  std::set<int> s;
  std::vector<bool> mark;
  std::vector<double> x; // Numerical values
  row_pattern(int n) { mark.resize(n, false); x.resize(n, 0); }
  bool empty() { return s.empty(); }
  int top() { return *s.begin(); }
  // Element access
  double& operator[](int i) { insert_pattern(i); return x[i]; }
  void erase(int i) { erase_pattern(i); x[i] = 0; }
  // Private
  void insert_pattern(int i) { if (!mark[i]) { s.insert(i); mark[i] = true; } }
  void erase_pattern (int i) { if (mark[i]) { s.erase(i); mark[i] = false; } }
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
   - Output is a lower triangular matrix representing pattern and numerical values of LDL fcatorization: 'L-I+D'
   - Factor entries below tolerance 'tol' are dropped *except* entries that are already present in A.
   - Passing 'tol=0' gives a complete factorization.
*/
cs* cs_ichol (const cs *C, double tol)
{
    double d, lki, *Lx, *Cx ;
    csi i, p, k, n, *Li, *Lp, *Cp, *Ci ;
    cs *L ;
    n = C->n ;
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    row_pattern x(n);
    std::vector<col_pattern > Lcol(n);
    std::vector<double> Ldiag(n);
    std::vector<bool> must_keep(n, false);
    for (k = 0 ; k < n ; k++)       /* compute L(k,:) for L*L' = C */
    {
        /* --- Nonzero pattern of L(k,:) ------------------------------------ */
        x.erase(k) ;                                /* x (0:k) is now zero */
        for (p = Cp [k] ; p < Cp [k+1] ; p++)       /* x = full(triu(C(:,k))) */
        {
            if (Ci [p] <= k) {
                x [Ci [p]] = Cx [p] ;
                must_keep[Ci[p]] = true;
            }
        }
        d = x [k] ;                     /* d = C(k,k) */
        x.erase(k);
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
        //if (d <= 0) return (cs_ndone (N, E, c, x_, 0)) ; /* not pos def */
        Ldiag[k] = d ;   /* store L(k,k) = sqrt (d) in column k */
    }
    /* finalize L */
    int nzmax=0;
    for (auto x : Lcol) nzmax += x.size() + 1;
    L = cs_spalloc (n, n, nzmax, 1, 0) ;    /* allocate result */
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
     of L.

     Note: 'cs_ichol_update' needs the row pattern (R = Transposed of
     L) without numerical values (R -> x). This can in principle be
     precomputed, but for now we skip this optimization and calculate
     R via cs_transpose inside the function.
*/
void cs_ichol_update (const cs *A, cs *L)
{
    double d, lki, *Lx, *x, *Cx ;
    csi top, i, p, k, n, *Li, *Lp, *Ri, *Rp, *c, *Cp, *Ci ;
    cs *C;
    n = A->n ;
    c = CS_MALLOC (2*n, SIZEOF (csi)) ;     /* get csi workspace */
    x = CS_MALLOC (n, SIZEOF (double)) ;    /* get double workspace */
    C = (cs*) A;
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    Lp = L->p ; Li = L->i ; Lx = L->x ;
    cs* R = cs_transpose(L, 0); // Pattern only
    Rp = R->p ; Ri = R->i ;
    for (k = 0 ; k < n ; k++) c [k] = Lp [k];
    double M = 0;
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
        /* solve L(0:k-1,0:k-1) * x = C(:,k) */
        for ( top=Rp[k]; top<Rp[k+1]-1; top++)
        {
            i = Ri [top] ;               /* s [top..n-1] is pattern of L(k,:) */
            double xi = x[i];
            x [i] = 0 ;                 /* clear x for k+1st iteration */
            lki = xi / Lx [Lp [i]] ;
            for (p = Lp [i] + 1 ; p < c [i] ; p++)
            {
                x [Li [p]] -= Lx [p] * xi ;
            }
            d -= lki * lki * Lx [Lp [i]] ;            /* d = d - L(k,i)*L(k,i) */
            p = c [i]++ ;
            Li [p] = k ;                /* store L(k,i) in column i */
            Lx [p] = lki ;
        }
        // Clear just the row pattern
        for ( top=Rp[k]; top<Rp[k+1]-1; top++) {
          x[Ri[top]] = 0;
        }
        // Now clear entire column union while recording the max
        for ( top=Rp[k]; top<Rp[k+1]-1; top++) {
          i = Ri [top] ;
          for (p = Lp [i] + 1 ; p < c [i] ; p++)
            {
              double tmp = x[Li [p]];
              if (tmp != 0) {
                M = std::max(M, std::abs(tmp / d));
              }
              x [Li [p]] = 0;
            }
        }
        /* --- Compute L(k,k) ----------------------------------------------- */
        //if (d <= 0) return (cs_ndone (N, E, c, x, 0)) ; /* not pos def */
        p = c [k]++ ;
        Li [p] = k ;                /* store L(k,k) = sqrt (d) in column k */
        Lx [p] = d ;
    }
    std::cout << "M = " << M << "\n";
    //Lp [n] = cp [n] ;               /* finalize L */
    cs_spfree(R);
    cs_free(c);
    cs_free(x);
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
  for (int col = nc; col > 0; ) {
    col--;
    // Get d
    int p_diag = Lp[col];
    double d = Lx[p_diag];
    // scatter this off-diagonal column
    for (int p = Lp[col] + 1; p < Lp[col+1]; p++) {
      w[Li[p]] = Lx[p];
    }
    // Loop through target subset (non-zeros of relevant columns)
    for (int p = Lp[col] + 1; p < Lp[col+1]; p++) {
      int tcol = Li[p]; // target column
      for (int tp = Lp[tcol]; tp < Lp[tcol+1]; tp++) {
        int trow = Li[tp]; // target row
        if (w[trow] != 0 && w[tcol] != 0) {
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
cs *cs_transpose (const cs *A, csi values)
{
    csi p, q, j, *Cp, *Ci, n, m, *Ap, *Ai, *w ;
    double *Cx, *Ax ;
    cs *C ;
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    C = CS_SPALLOC (n, m, Ap [n], values && Ax, 0) ;       /* allocate result */
    w = CS_CALLOC (m, SIZEOF (csi)) ;                      /* get workspace */
    if (!C || !w) return (cs_done (C, w, NULL, 0)) ;       /* out of memory */
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (p = 0 ; p < Ap [n] ; p++) w [Ai [p]]++ ;          /* row counts */
    cs_cumsum (Cp, w, m) ;                                 /* row pointers */
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            Ci [q = w [Ai [p]]++] = j ; /* place A(i,j) as entry C(j,i) */
            if (Cx) Cx [q] = Ax [p] ;
        }
    }
    return (cs_done (C, w, NULL, 1)) ;  /* success; free w and return C */
}
/* allocate a sparse matrix (triplet form or compressed-column form) */
cs *cs_spalloc (csi m, csi n, csi nzmax, csi values, csi triplet)
{
    cs *A = CS_CALLOC (1, SIZEOF (cs)) ;    /* allocate the cs struct */
    if (!A) return (NULL) ;                 /* out of memory */
    A->m = m ;                              /* define dimensions and nzmax */
    A->n = n ;
    A->nzmax = nzmax = CS_MAX (nzmax, 1) ;
    A->nz = triplet ? 0 : -1 ;              /* allocate triplet or comp.col */
    A->p = CS_MALLOC (triplet ? nzmax : n+1, SIZEOF (csi)) ;
    A->i = CS_MALLOC (nzmax, SIZEOF (csi)) ;
    A->x = values ? CS_MALLOC (nzmax, SIZEOF (double)) : NULL ;
    return ((!A->p || !A->i || (values && !A->x)) ? cs_spfree (A) : A) ;
}

/* free a sparse matrix */
cs *cs_spfree (cs *A)
{
    if (!A) return (NULL) ;     /* do nothing if A already NULL */
    cs_free (A->p) ;
    cs_free (A->i) ;
    cs_free (A->x) ;
    return ((cs *) cs_free (A)) ;   /* free the cs struct and return NULL */
}

/* free workspace and return a sparse matrix result */
cs *cs_done (cs *C, void *w, void *x, csi ok)
{
    cs_free (w) ;                       /* free workspace */
    cs_free (x) ;
    return (ok ? C : cs_spfree (C)) ;   /* return result if OK, else free it */
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
  int* Dim = INTEGER(ALLOC_SLOT(ans, install("Dim"), INTSXP, 2));
  Dim[0] = L->m; Dim[1] = L->n;
  /* p : */
  int* p = INTEGER(ALLOC_SLOT(ans, install("p"), INTSXP, L->n + 1));
  Memcpy(p, (int*)(L->p), L->n + 1);
  /* i : */
  int* i = INTEGER(ALLOC_SLOT(ans, install("i"), INTSXP, L->nzmax));
  Memcpy(i, (int*)(L->i), L->nzmax);
  /* x : */
  if (L->x != NULL) {
    double* x = REAL(ALLOC_SLOT(ans, install("x"), REALSXP, L->nzmax));
    Memcpy(x, (int*)(L->x), L->nzmax);
  }
}
cs r2cs(SEXP X) {
  cs A;
  A.nzmax = LENGTH(R_do_slot(X, install("x")));
  A.m = INTEGER(R_do_slot(X, install("Dim")))[0];
  A.n = INTEGER(R_do_slot(X, install("Dim")))[1];
  A.p = INTEGER(R_do_slot(X, install("p")));
  A.i = INTEGER(R_do_slot(X, install("i")));
  bool has_x = R_has_slot(X, install("x"));
  A.x = has_x ? REAL(R_do_slot(X, install("x"))) : NULL;
  A.nz = -1;
  return A;
}
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
  cs* L = cs_ichol(&A, chol_tol);
  SEXP claes = PROTECT(R_do_MAKE_CLASS("dgCMatrix"));  // OR dgTMatrix ?
  SEXP ans = PROTECT(R_do_new_object(claes));
  cs2r(ans, L);
  cs_spfree(L);
  UNPROTECT(2);
  return ans;
}
// Update incomplete Cholesky factor L
extern "C"
SEXP tmb_ichol_update(SEXP X, SEXP Y) {
  cs A = r2cs(X);
  cs L = r2cs(Y);
  cs_ichol_update(&A, &L);
  return Y;
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
