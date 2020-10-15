#ifdef _OPENMP

/* Override R-API with thread safe versions.

   FIXME: Still some missing e.g. XLENGTH

   FIXME: To minimize overhead one should use as few R-API calls as
   possible, i.e. avoid doing REAL(x)[i] in a loop.
*/

inline SEXP Ts_getAttrib(SEXP x, SEXP y) {
  SEXP ans;
#pragma omp critical
  {
    ans = Rf_getAttrib(x, y);
  }
  return ans;
}

inline SEXP Ts_STRING_ELT(SEXP x, size_t i) {
  SEXP ans;
#pragma omp critical
  {
    ans = STRING_ELT(x, i);
  }
  return ans;
}

inline const char* Ts_CHAR(SEXP x) {
  const char* ans;
#pragma omp critical
  {
    ans = R_CHAR(x);
  }
  return ans;
}

inline SEXP Ts_VECTOR_ELT(SEXP x, size_t i) {
  SEXP ans;
#pragma omp critical
  {
    ans = VECTOR_ELT(x, i);
  }
  return ans;
}

inline R_len_t Ts_length(SEXP x) {
  R_len_t ans;
#pragma omp critical
  {
    ans = Rf_length(x);
  }
  return ans;
}

inline int* Ts_INTEGER(SEXP x) {
  int* ans;
#pragma omp critical
  {
    ans = INTEGER(x);
  }
  return ans;
}

inline double* Ts_REAL(SEXP x) {
  double* ans;
#pragma omp critical
  {
    ans = REAL(x);
  }
  return ans;
}

extern "C"
inline void Ts_GetRNGstate() {
#pragma omp master
  {
    GetRNGstate();
  }
  // Threads wait for master
#pragma omp barrier
}

inline Rboolean Ts_isNumeric(SEXP x) {
  Rboolean ans;
#pragma omp critical
  {
    ans = Rf_isNumeric(x);
  }
  return ans;
}

inline int Ts_LENGTH(SEXP x) {
  int ans;
#pragma omp critical
  {
    ans = LENGTH(x);
  }
  return ans;
}

inline SEXP Ts_install(const char *x) {
  SEXP ans;
#pragma omp critical
  {
    ans = Rf_install(x);
  }
  return ans;
}

// Redefine
#define Rf_getAttrib   Ts_getAttrib
#define STRING_ELT     Ts_STRING_ELT
#undef  CHAR
#define CHAR           Ts_CHAR
#define VECTOR_ELT     Ts_VECTOR_ELT
#define Rf_length      Ts_length
#define INTEGER        Ts_INTEGER
#define REAL           Ts_REAL
#define GetRNGstate    Ts_GetRNGstate
#define Rf_isNumeric   Ts_isNumeric
#define LENGTH         Ts_LENGTH
#define Rf_install     Ts_install

#endif // _OPENMP
