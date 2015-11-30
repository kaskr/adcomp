// Copyright (C) 2013-2015 Kasper Kristensen
// License: GPL-2

# include <R.h>
# include <Rinternals.h>
# ifdef _OPENMP
#include <omp.h>
# endif

/* openmp controller */
SEXP omp_num_threads(SEXP x) {
#ifdef SUPPORT_OPENMP
  if(!isNull(x)){
    int n=INTEGER(x)[0];
    omp_set_num_threads(n);
  }
  SEXP ans;
  PROTECT(ans = allocVector(INTSXP,1));
  INTEGER(ans)[0]=omp_get_max_threads();
  UNPROTECT(1);
  return ans;
#else
  error("Openmp not supported.");
#endif
}

/* Avoid S4 overhead when changing x-slot:
   Set xslot to SEXP pointer i.e. x@x <- y
 */
SEXP setxslot(SEXP x, SEXP y){
  setAttrib(x,install("x"),y);
  return x;
}

/* Is external pointer nil ? */
SEXP isNullPointer(SEXP pointer) {
  return ScalarLogical(!R_ExternalPtrAddr(pointer));
}

/* Mark factor with numerical values of hessian */
SEXP tmb_mark_factor(SEXP L, SEXP H){
  setAttrib(L,
	    install("Hx"),
	    duplicate(getAttrib(H, install("x"))));
  return R_NilValue;
}

/* Check if L already contains the factorization of H */
SEXP tmb_match_factor(SEXP L, SEXP H){
  SEXP Hx = getAttrib(H, install("x"));
  SEXP Lx = getAttrib(L, install("Hx"));
  if(Lx == R_NilValue)return ScalarLogical(0);
  int n = LENGTH(Lx);
  if(n != LENGTH(Hx)) error("'tmb_match_factor' unequal lengths");
  for(int i=0; i < n; i++){
    if (REAL(Lx)[i] != REAL(Hx)[i])
      return ScalarLogical(0);
  }
  return ScalarLogical(1);
}
