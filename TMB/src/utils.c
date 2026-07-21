// Copyright (C) 2013-2015 Kasper Kristensen
// License: GPL-2

# include <R.h>
# include <Rinternals.h>

/* Avoid S4 overhead when changing x-slot:
   Set xslot to SEXP pointer i.e. x@x <- y
 */
SEXP setxslot(SEXP x, SEXP y){
  setAttrib(x,install("x"),y);
  return x;
}
/* Set slots inplace */
SEXP setslot(SEXP x, SEXP nm, SEXP y) {
  SEXP symb = Rf_installChar(STRING_ELT(nm, 0));
  Rf_setAttrib(x, symb, y);
  return x;
}

/* Is external pointer nil ? */
SEXP isNullPointer(SEXP pointer) {
  return ScalarLogical(TYPEOF(pointer) == EXTPTRSXP && !R_ExternalPtrAddr(pointer));
}
