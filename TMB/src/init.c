#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "Matrix.h"
cholmod_common c;

#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}

SEXP omp_num_threads(SEXP x);
SEXP isNullPointer(SEXP pointer);
SEXP setxslot(SEXP x, SEXP y);
SEXP tmb_invQ(SEXP Lfac);
SEXP tmb_invQ_tril_halfdiag(SEXP Lfac);
SEXP match_pattern(SEXP A_, SEXP B_);
SEXP tmb_sparse_izamd(SEXP A_, SEXP mark_, SEXP diag_);
SEXP tmb_half_diag(SEXP A_);
SEXP have_tmb_symbolic(void);
SEXP tmb_symbolic(SEXP Qp);
SEXP tmb_destructive_CHM_update(SEXP L, SEXP H, SEXP mult);
SEXP tmb_CHMfactor_solve(SEXP L_, SEXP y_);

static R_CallMethodDef CallEntries[] = {
    CALLDEF(omp_num_threads, 1),
    CALLDEF(isNullPointer, 1),
    CALLDEF(setxslot, 2),
    CALLDEF(tmb_invQ, 1),
    CALLDEF(tmb_invQ_tril_halfdiag, 1),
    CALLDEF(match_pattern, 2),
    CALLDEF(tmb_sparse_izamd, 3),
    CALLDEF(tmb_half_diag, 1),
    CALLDEF(have_tmb_symbolic, 0),
    CALLDEF(tmb_symbolic, 1),
    CALLDEF(tmb_destructive_CHM_update, 3),
    CALLDEF(tmb_CHMfactor_solve, 2),
    {NULL, NULL, 0}
};

void R_init_TMB(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, (Rboolean)FALSE);
    M_R_cholmod_start(&c);
    c.error_handler = NULL; // Disable CHOLMOD warnings
}

void R_unload_TMB(DllInfo *dll)
{
    M_cholmod_finish(&c);
}
