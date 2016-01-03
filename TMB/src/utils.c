// Copyright (C) 2013-2015 Kasper Kristensen
// License: GPL-2

# include <R.h>
# include <Rinternals.h>
# ifdef _OPENMP
#include <omp.h>
# endif

/* openmp controller */
SEXP omp_num_threads(SEXP x) {
#ifdef _OPENMP
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
