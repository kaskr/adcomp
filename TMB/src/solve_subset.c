// Copyright (C) 2013-2015 Kasper Kristensen
// License: GPL-2

/* ========================================================== 
   Supernodal version of solvesubset for CHOLMOD supernodal
   sparse Cholesky structures.

   Description:
   * Given the factorization A=LL'.
   * Calculate the inverse S=A^-1 on the sparseness pattern
     of LL'.
   NOTE: In the dense case this is equivalent with "DPOTRI".

   Algorithm (Recursions):
   * s = indices of supernode
   * p = non-zero indices below supernode

   [ L(s,s) L(s,p) ]        [ S(s,s) S(s,p) ]
   [ L(p,s) L(p,p) ]        [ S(p,s) S(p,p) ]

   1. S(s,p) = -L(s,s)^t^-1 * L(s,p) * S(p,p)
   2. S(s,s) = -L(s,s)^t^-1 * L(s,p) * S(p,s) + (L(s,s)*L(s,s)^t)^-1

   Rewritten:
   M0 = (L(s,s)*L(s,s)^t)^-1 (DPOTRI)
   M = -L(p,s) * L(s,s)^-1 (DTRSM)
   S(p,s) =  S(p,p) * M (DSYMM)
   S(s,s) =  M^t * S(p,s) + M0 (DGEMM)

   ====> IMPLEMENTATION:
   Compute dense submatrix Lss=L(s,s)
   If p not empty {
     1. Compute M = -L(p,s) * L(s,s)^-1 :
        Compute dense submatrix M=L(p,s)
        M = - M * Lss^-1
        DTRSM("R", "L",np,ns,-1.0,Lss,ns,M,np);
     2. Compute S(p,s) =  S(p,p) * M
        Compute dense submatrix Spp=S(p,p)
     3. Compute S(s,s) =  M^t * S(p,s) + M0
   }

   M0 = (L(s,s)*L(s,s)^t)^-1 = DPOTRI("U",dims,Lss,dims,&info);
  
   ========================================================== 
*/

#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "Matrix.h"

/* Copy-pasted from "Writing R Extensions". A similar spell is present
   in Matrix.h so might not be needed anymore ? */
#ifdef __GNUC__
// this covers gcc, clang, icc
# undef alloca
# define alloca(x) __builtin_alloca((x))
#elif defined(HAVE_ALLOCA_H)
// needed for native compilers on Solaris and AIX
# include <alloca.h>
#endif

extern cholmod_common c; // See init.c

/* Temporary fix until 'AS_CHM_FR__' is part of Matrix API.
   In future replace code below by something like:

#ifdef AS_CHM_FR__
#undef AS_CHM_FR
#define AS_CHM_FR AS_CHM_FR__
#endif

*/
#undef AS_CHM_FR
#define AS_CHM_FR(x) tmb_as_cholmod_factor3((CHM_FR)alloca(sizeof(cholmod_factor)), x, FALSE)
CHM_FR tmb_as_cholmod_factor3(CHM_FR ans, SEXP x, Rboolean do_check)
{
  int *type = INTEGER(GET_SLOT(x, install("type")));
  SEXP tmp;
  memset(ans, 0, sizeof(cholmod_factor)); /* zero the struct */
  ans->itype = CHOLMOD_INT;	/* characteristics of the system */
  ans->dtype = CHOLMOD_DOUBLE;
  ans->z = (void *) NULL;
  ans->xtype = CHOLMOD_REAL;
  ans->ordering = type[0];	/* unravel the type */
  ans->is_ll = (type[1] ? 1 : 0);
  ans->is_super = (type[2] ? 1 : 0);
  ans->is_monotonic = (type[3] ? 1 : 0);
  SEXP Matrix_permSym = install("perm");
  tmp = GET_SLOT(x, Matrix_permSym);
  ans->minor = ans->n = LENGTH(tmp); ans->Perm = INTEGER(tmp);
  ans->ColCount = INTEGER(GET_SLOT(x, install("colcount")));
  ans->z = ans->x = (void *) NULL;
  SEXP Matrix_xSym = install("x");
  tmp = GET_SLOT(x, Matrix_xSym);
  ans->x = REAL(tmp);
  if (ans->is_super) {	/* supernodal factorization */
    ans->xsize = LENGTH(tmp);
    ans->maxcsize = type[4]; ans->maxesize = type[5];
    ans->i = (int*)NULL;
    tmp = GET_SLOT(x, install("super"));
    ans->nsuper = LENGTH(tmp) - 1; ans->super = INTEGER(tmp);
    tmp = GET_SLOT(x, install("pi"));
    ans->pi = INTEGER(tmp);
    tmp = GET_SLOT(x, install("px"));
    ans->px = INTEGER(tmp);
    tmp = GET_SLOT(x, install("s"));
    ans->ssize = LENGTH(tmp); ans->s = INTEGER(tmp);
  } else {
    Rf_error("Unexpected");
  }
  return ans;
}

SEXP tmb_destructive_CHM_update(SEXP L, SEXP H, SEXP mult)
{
    CHM_FR f = AS_CHM_FR(L);
    CHM_SP A = AS_CHM_SP__(H);
    double mm[2] = {0, 0};
    mm[0] = asReal(mult);
    // NB: Result depends if A is "dsC" or "dgC"; the latter case assumes we mean AA' !!!
    /*
      cholmod_factorize_p return value:

      TRUE:  CHOLMOD_OK, CHOLMOD_NOT_POSDEF, CHOLMOD_DSMALL

      FALSE: CHOLMOD_NOT_INSTALLED, CHOLMOD_OUT_OF_MEMORY,
             CHOLMOD_TOO_LARGE, CHOLMOD_INVALID, CHOLMOD_GPU_PROBLEM
    */
    if (!M_cholmod_factorize_p(A, mm, (int*)NULL, 0 /*fsize*/, f, &c))
      /* -> ./CHOLMOD/Cholesky/cholmod_factorize.c */
      error("cholmod_factorize_p failed: status %d, minor %d of ncol %d",
            c.status, f->minor, f->n);
    int ok = (c.status == CHOLMOD_OK);
    return ScalarLogical(ok);
}

SEXP tmb_CHMfactor_solve(SEXP L_, SEXP y_)
{
    CHM_FR L = AS_CHM_FR(L_);
    int n = LENGTH(y_);
    CHM_DN y = N_AS_CHM_DN(REAL(y_), n, 1);
    SEXP x = PROTECT( NEW_NUMERIC( n ) );
    CHM_DN sol = M_cholmod_solve(CHOLMOD_A, L, y, &c);
    memcpy(REAL(x), sol->x, n * sizeof(double));
    M_cholmod_free_dense(&sol, &c);
    UNPROTECT(1);
    return x;
}

// Notes about the CHOLMOD super-nodal storage. 
// According to the documentation of CHOLMOD we have:
// ==================================================
//               k1 = Super [s] ;            /* s contains columns k1 to k2-1 of L */
//               k2 = Super [s+1] ;
//               nscol = k2 - k1 ;           /* # of columns in all of s */
//               psi = Lpi [s] ;             /* pointer to first row of s in Ls */
//               psx = Lpx [s] ;             /* pointer to first row of s in Lx */
//               psend = Lpi [s+1] ;         /* pointer just past last row of s in Ls */
//               nsrow = psend - psi ;       /* # of rows in all of s */
//        * See discussion in cholmod_change_factor.c:
//        * (3) supernodal symbolic:  A representation of the nonzero pattern of the
//        *      supernodes for a supernodal factorization.  There are L->nsuper
//        *      supernodes.  Columns L->super [k] to L->super [k+1]-1 are in the kth
//        *      supernode.  The row indices for the kth supernode are in
//        *      L->s [L->pi [k] ... L->pi [k+1]-1].  The numerical values are not
//        *      allocated (L->x), but when they are they will be located in
//        *      L->x [L->px [k] ... L->px [k+1]-1], and the L->px array is defined
//        *      in this factor type.

					    
/* Extract dense block x[p,q] of sparse matrix x */
CHM_DN densesubmatrix(CHM_SP x, int *p, int np, int *q, int nq, cholmod_common *c){
  CHM_DN ans = M_cholmod_allocate_dense(np,nq,np,CHOLMOD_REAL,c);
  double *w = malloc(x->nrow*sizeof(double));
  int *xi=x->i;
  int *xp=x->p;
  double *xx=x->x;
  double *ansx=ans->x;
  int col, row;
  for(int j=0;j<nq;j++){
    col=q[j];
    /* scatter x[:,col] over w */
    for(int i=xp[col];i<xp[col+1];i++){
      w[xi[i]]=xx[i];
    }
    /* Copy w[p] to ans[:,col] */
    for(int i=j;i<np;i++){
      row=p[i];
      ansx[i+j*np]=w[row];
    }
  }
  free(w);
  return ans;
}

/* Perform recursions for k'th supernode */
void tmb_recursion_super(CHM_SP Lsparse, int k, CHM_FR L, cholmod_common *c){
  int* super=L->super;
  int* Ls=L->s;
  int* Lpi=L->pi;
  int ncol=super[k+1]-super[k]; /* ncol of supernode */
  int nrow=Lpi[k+1]-Lpi[k]; /* Number of rows in supernode */
  /* q contains row-indices of *entire* supernode */
  /* p contains row-indices excluding those of diagonal */
  /* s contains row-indices of diagonal - setdiff(q,p) */
  int* q=Ls+Lpi[k];    /* Pointer to L->s [L->pi [k]] */
  int nq=nrow;         /* length of q */
  // int* p=q+ncol;       /* Exclude triangle in diagonal */
  int np=nq-ncol;      /* length of p */
  int* s=q;
  int ns=ncol;         /* length of s */
  /* do not sort because p is sorted */
  int info; /* For lapack */
  int i,j;
  double ONE=1.0, ZERO=0.0, MONE=-1.0;
  CHM_DN x = densesubmatrix(Lsparse,q,nq,q,nq,c);
  double *xx=x->x;
  double *Lss=xx, *Lps=xx+ns, *Ssp=xx+(nq*ns), *Spp=xx+(nq*ns+ns);
  /* Workspace to hold output from dsymm */
  double *wrk=malloc(nq*ns*sizeof(double));
  double *wrkps=wrk+ns;
  if(np>0){
    F77_CALL(dtrsm)("Right","Lower","No transpose","Not unit",
		    &np,&ns,&MONE,Lss,&nq,Lps,&nq);
    for(i=ns;i<nq;i++){for(j=0;j<ns;j++)Lss[j+nq*i] = Lss[i+nq*j];} /* Copy Transpose*/
    memcpy(wrk,Lss,nq*ns*sizeof(double));    
    F77_CALL(dsymm)("Left","Lower",&np,&ns,&ONE,Spp,&nq,Lps,&nq,&ZERO,wrkps,&nq);
    memcpy(Lss,wrk,nq*ns*sizeof(double));
    F77_CALL(dpotri)("L",&ns,Lss,&nq,&info);
    F77_CALL(dgemm)("N","N",&ns,&ns,&np,&ONE,Ssp,&nq,Lps,&nq,&ONE,Lss,&nq);
  } else {
    F77_CALL(dpotri)("L",&ns,Lss,&nq,&info);
  }
  /* ----------- Fill results into L(q,s) -------*/
  double *Lx=Lsparse->x;
  int *Lp=Lsparse->p;
  int m=Lp[s[0]];
  for(j=0;j<ns;j++){
    for(i=j;i<nq;i++){
      Lx[m]=Lss[i+j*nq];
      m++;
    }
  }

  /* Count flops */
  /* if(control.countflops){ */
  /*   double Ns=ns; */
  /*   double Np=np; */
  /*   if(np>0){ */
  /*     flopcount[0]+=(Ns*(Ns+1))*0.5*Np; /\* dtrsm *\/ */
  /*     flopcount[1]+=Np*Np*Ns; /\* dsymm *\/ */
  /*     flopcount[2]+=2.0*(Ns*Ns*Ns)/3.0; /\* dpotri *\/ */
  /*     flopcount[3]+=Ns*Np*Ns; /\* dgemm *\/ */
  /*   } else { */
  /*     flopcount[2]+=2.0*(Ns*Ns*Ns)/3.0; /\* dpotri *\/ */
  /*   } */
  /* } */

  /* Clean up */
  M_cholmod_free_dense(&x,c);
  free(wrk);
}

CHM_SP tmb_inv_super(CHM_FR Lfac, cholmod_common *c){

  /* Convert factor to sparse without modifying factor */
  CHM_FR Ltmp = M_cholmod_copy_factor(Lfac,c);
  CHM_SP L = M_cholmod_factor_to_sparse(Ltmp,c);
  M_cholmod_free_factor(&Ltmp,c);

  /* Loop over supernodes */
  int nsuper=Lfac->nsuper;
  for(int k=nsuper-1;k>=0;k--)tmb_recursion_super(L,k,Lfac,c);

  /* Change to symm lower */
  L->stype=-1; 
  return L;
}

SEXP tmb_invQ(SEXP Lfac){
  CHM_FR L=AS_CHM_FR(Lfac);
  CHM_SP iQ = tmb_inv_super(L, &c);
  return M_chm_sparse_to_SEXP(iQ, 1 /* Free */ , 0, 0, "", R_NilValue);
}

void half_diag(CHM_SP A){
  int ncol=A->ncol;
  double *Ax;
  int *Ai, *Ap, i;
  Ai=A->i; Ap=A->p; Ax=A->x;
  for(int j=0;j<ncol;j++){
    for(int k=Ap[j];k<Ap[j+1];k++){
      i=Ai[k];
      if(i==j)Ax[k]=.5*Ax[k];
    }
  }
}

SEXP tmb_invQ_tril_halfdiag(SEXP Lfac){
  CHM_FR L=AS_CHM_FR(Lfac);
  CHM_SP iQ = tmb_inv_super(L, &c);
  half_diag(iQ);
  iQ->stype=0; /* Change to non-sym */
  return M_chm_sparse_to_SEXP(iQ, 1 /* Free */ , -1 /* uplo="L" */ , 0, "", R_NilValue);
}

/* Given sparse matrices A and B (sorted columns).
   Assume pattern of A is a subset of pattern of B.
   (This also includes cases where dimension of B larger than dim of A)
   Return integer vector p of same length as A@x such that
     " A@i == B@i[p] and A@j == B@j[p] "
*/
SEXP match_pattern(SEXP A_, SEXP B_){
  CHM_SP A=AS_CHM_SP(A_);
  CHM_SP B=AS_CHM_SP(B_);
  int *Ai=A->i, *Bi=B->i, *Ap=A->p, *Bp=B->p;
  int ncol=A->ncol,i,j,k;
  int index; // index match
  SEXP ans;
  if(A->ncol>B->ncol)error("Must have dim(A)<=dim(B)");
  PROTECT(ans=NEW_INTEGER(A->nzmax));
  int *pans=INTEGER(ans);
  for(j=0;j<ncol;j++){
    index=Bp[j]; // Start at top of B(:,j)
    for(k=Ap[j];k<Ap[j+1];k++){
      i=Ai[k];
      for(;Bi[index]!=i;index++){ // Find next match
	if(index>=Bp[j+1]){
	  UNPROTECT(1);
	  error("No match");
	}
      }
      *pans=index+1; pans++; // R-index !
    }  
  }
  UNPROTECT(1);
  return ans;
}

/*
  Sparse version of 'insert zeros and modify diagonal' ( izamd )
  keeping the sparsity pattern un-modified:

  A[:,p] <- 0;  A[p,:] <- 0;  diag(A)[p] <- d;

  A_      : Sparse matrix to modify
  mark_   : Logical (int) index vector of rows and columns.
  diag_   : Diagonal replacement (double).

*/
SEXP tmb_sparse_izamd(SEXP A_, SEXP mark_, SEXP diag_){
  CHM_SP A = AS_CHM_SP(A_);
  int *Ai=A->i, *Ap=A->p;
  double *Ax = A->x;
  int ncol=A->ncol;
  int *mark = INTEGER(mark_);
  double diag = REAL(diag_)[0];
  int i, l=0;
  for(int j = 0; j < ncol; j++){
    for(int k = Ap[j]; k < Ap[j+1]; k++){
      i = Ai[k];
      if (mark[i]) Ax[l] = 0;
      if (mark[j]) Ax[l] = 0;
      if ( (mark[i] || mark[j]) && (i == j) ) Ax[l] = diag;
      l++;
    }
  }
  return A_;
}

/* Half the diagonal of a matrix (note: modifies input) */
SEXP tmb_half_diag(SEXP A_){
  CHM_SP A = AS_CHM_SP(A_);
  half_diag(A);
  return A_;
}
