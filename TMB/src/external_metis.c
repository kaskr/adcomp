// Copyright (C) 2013-2015 Kasper Kristensen
// License: GPL-2

#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "Matrix.h"
#include <Rdefines.h>

#ifdef _USE_EXTERNAL_CHOLMOD_LIB_
/* ========================================================================== */
/* Error handler */
/* ========================================================================== */
void tmb_cholmod_error(int status, const char *file, int line, const char *message)
{
  /* From CHOLMOD/Include/cholmod_core.h : ...status values.
     zero means success, negative means a fatal error, positive is a warning.
  */
  warning(("Cholmod warning '%s' at file:%s, line %d"),
	  message, file, line);
}

cholmod_factor *cholmod_analyze
(
    /* ---- input ---- */
    cholmod_sparse *A, /* matrix to order and analyze */
    /* --------------- */
    cholmod_common *Common
) ;
int cholmod_free_factor
(
    /* ---- in/out --- */
    cholmod_factor **LHandle,	/* factor to free, NULL on output */
    /* --------------- */
    cholmod_common *Common
) ;

/* ========================================================================== */
/* Run the symbolic analysis and prepare workspace - only run once !!! */
/* ========================================================================== */
SEXP tmb_symbolic(SEXP Qp){
  cholmod_common c;
  M_R_cholmod_start(&c);
  /* TODO: More control from R */
  c.nmethods=9; 
  c.supernodal = CHOLMOD_SUPERNODAL;
  c.final_ll = TRUE;
  /*  Return quickly if not positive definite */
  c.quick_return_if_not_posdef=TRUE;
  c.error_handler=tmb_cholmod_error;
  int trace=1;

  CHM_SP Q = M_cholmod_copy(AS_CHM_SP(Qp), -1 /* symmetric lower */, 1 /*values*/, &c);
  CHM_FR LQ;
  CHM_FR LQ2;

  // Step 1: Run symbolic analysis with external cholmod library:
  if(trace)Rprintf("Entering externallib \n");
  c.itype=CHOLMOD_INT;
  LQ = cholmod_analyze(Q, &c); /* get fill-reducing permutation */
  if(trace)Rprintf("cholmod_analyze: status=%d \n",c.status);
  if(trace)Rprintf("Chosen ordering %d \n", c.selected);
  
  // Step 2: Grab the permutation:
  int *perm=LQ->Perm;

  // Step 3: Run symbolic analysis again, now with known permutation
  //         using the R cholmod interface routines
  if(trace)Rprintf("Running symbolic analysis \n");
  if(trace)Rprintf("User permutation \n");
  c.nmethods=1; 
  LQ2 = M_cholmod_analyze_p(Q, perm, NULL, 0, &c);
  cholmod_free_factor(&LQ,&c); // LQ Not needed anymore
  
  if(trace)Rprintf("Chosen ordering %d \n", c.selected);
  if(trace)Rprintf("Length of supernodal xslot %d \n", LQ2->xsize);
  if(trace)Rprintf("Flopcount %f \n", c.fl);
  double nnzL = LQ2->xsize;
  double nnzQ = Q->nzmax;
  if(trace)Rprintf("Fill-in ratio (nnz(L)/nnz(Q)) %f \n", nnzL/nnzQ);
  if(trace)Rprintf("Factor xtype %d \n", LQ2->xtype);

  // Step 4: Make sure factor has numerical values 
  if(trace)Rprintf("Running numerical factorization \n");
  M_cholmod_factorize(Q, LQ2, &c);
  if(trace)Rprintf("Done \n");
  
  // Cleanup
  M_cholmod_free_sparse(&Q,&c);
  M_cholmod_finish(&c);

  return M_chm_factor_to_SEXP(LQ2, 1 /* Free */);
}

#else

SEXP tmb_symbolic(SEXP Qp) {
  return R_NilValue;
}

#endif

SEXP have_tmb_symbolic(void) {
  SEXP ans;
  PROTECT(ans = NEW_INTEGER(1));
#ifdef _USE_EXTERNAL_CHOLMOD_LIB_
  INTEGER(ans)[0]=1;
#else
  INTEGER(ans)[0]=0;
#endif
  UNPROTECT(1);
  return ans;
}
