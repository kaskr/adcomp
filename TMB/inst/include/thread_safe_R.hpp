// getAttrib
// STRING_ELT
// CHAR
// VECTOR_ELT
// Rf_length
// INTEGER
// REAL
// GetRNGstate

SEXP Ts_getAttrib(SEXP x, SEXP y) {
  SEXP ans;
#pragma omp critical
  {
    ans = Rf_getAttrib(x, y);
  }
  return ans;
}

void Ts_GetRNGstate() {
#pragma omp critical
  {
    GetRNGstate();
  }
}
