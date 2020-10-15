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

SEXP Ts_STRING_ELT(SEXP x, size_t i) {
  SEXP ans;
#pragma omp critical
  {
    ans = STRING_ELT(x, i);
  }
  return ans;
}

void Ts_GetRNGstate() {
#pragma omp critical
  {
    GetRNGstate();
  }
}
