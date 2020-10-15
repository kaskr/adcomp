// Rf_getAttrib
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

const char* Ts_CHAR(SEXP x) {
  const char* ans;
#pragma omp critical
  {
    ans = R_CHAR(x);
  }
  return ans;
}

SEXP Ts_VECTOR_ELT(SEXP x, size_t i) {
  SEXP ans;
#pragma omp critical
  {
    ans = VECTOR_ELT(x, i);
  }
  return ans;
}

extern "C"
void Ts_GetRNGstate() {
#pragma omp critical
  {
    GetRNGstate();
  }
}

// Redefine
#define Rf_getAttrib   Ts_getAttrib
#define STRING_ELT     Ts_STRING_ELT
#undef  CHAR
#define CHAR           Ts_CHAR
#define VECTOR_ELT     Ts_VECTOR_ELT
#define GetRNGstate    Ts_GetRNGstate
