// Copyright (C) 2013-2015 Kasper Kristensen
// License: GPL-2

# ifdef _OPENMP
#include <omp.h>
# endif
# include <R.h>
# include <Rinternals.h>

/* openmp controller */
SEXP omp_num_threads(SEXP x) {
#ifdef _OPENMP
  if( !isNull(x) ){
    int n = INTEGER(x)[0];
    omp_set_num_threads( n );
  }
  return ScalarInteger( omp_get_max_threads() );
#else
  warning("OpenMP not supported.");
  return ScalarInteger( 0 );
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
