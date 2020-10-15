// getAttrib
// STRING_ELT
// CHAR
// VECTOR_ELT
// Rf_length
// INTEGER
// REAL
// GetRNGstate

Ts_getAttrib(SEXP x, SEXP y) {
#pragma omp critical
  {
    SEXP ans = Rf_getAttrib(x, y);
  }
  return ans;
}

void Ts_GetRNGstate() {
#pragma omp critical
  {
    GetRNGstate();
  }
}
