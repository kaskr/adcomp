#ifndef HAVE_CONFIG_HPP
#define HAVE_CONFIG_HPP

/* ========================================================================== */
/* === Configuration ======================================================== */
/* ========================================================================== */

// Begin Configuration
#ifndef TMBAD_SCALAR_TYPE
#define TMBAD_SCALAR_TYPE double
#endif
#ifndef TMBAD_INDEX_TYPE
#define TMBAD_INDEX_TYPE unsigned int
#endif
#ifndef TMBAD_INDEX_VECTOR
#define TMBAD_INDEX_VECTOR std::vector<TMBAD_INDEX_TYPE>
#endif
#ifndef TMBAD_REPLAY_TYPE
#define TMBAD_REPLAY_TYPE ad_aug
#endif
#ifndef TMBAD_MAX_NUM_THREADS
#define TMBAD_MAX_NUM_THREADS 48
#endif
#ifndef TMBAD_COMPRESS_TOL
#define TMBAD_COMPRESS_TOL 16
#endif
#ifndef TMBAD_HASH_TYPE
#define TMBAD_HASH_TYPE unsigned int
#endif
#ifndef TMBAD_DETERMINISTIC_HASH
#define TMBAD_DETERMINISTIC_HASH true
#endif
#ifndef TMBAD_UNION_OR_STRUCT
#define TMBAD_UNION_OR_STRUCT union
#endif
#ifndef TMBAD_MIN_PERIOD_REP
#define TMBAD_MIN_PERIOD_REP 10
#endif
#ifdef __cpp_constexpr
#define CONSTEXPR constexpr
#else
#define CONSTEXPR
#endif
// End Configuration

#ifdef _OPENMP
#include <omp.h>
#define TMBAD_THREAD_NUM omp_get_thread_num()
#define TMBAD_SHARED_PTR TMBad::omp_shared_ptr
#else
#define TMBAD_SHARED_PTR std::shared_ptr
#define TMBAD_THREAD_NUM 0
#endif

/* ========================================================================== */
/* === Common macros ======================================================== */
/* ========================================================================== */

#ifndef ASSERT
#define ASSERT(x)                                \
  if (!(x)) {                                    \
    Rcerr << "ASSERTION FAILED: " << #x << "\n"; \
    abort();                                     \
  }
#define ASSERT2(x, msg)                          \
  if (!(x)) {                                    \
    Rcerr << "ASSERTION FAILED: " << #x << "\n"; \
    Rcerr << "POSSIBLE REASON: " << msg << "\n"; \
    abort();                                     \
  }
#endif
/* Test whether the maximum of TMBAD_INDEX_TYPE has been exceeded */
#define TMBAD_INDEX_OVERFLOW(x) \
  ((size_t)(x) >= (size_t)std::numeric_limits<TMBAD_INDEX_TYPE>::max())
#define xstringify(s) stringify(s)
#define stringify(s) #s

// Avoid typing nightmare.
#define INHERIT_CTOR(A, B)                                       \
  A() {}                                                         \
  template <class T1>                                            \
  A(const T1 &x1) : B(x1) {}                                     \
  template <class T1, class T2>                                  \
  A(const T1 &x1, const T2 &x2) : B(x1, x2) {}                   \
  template <class T1, class T2, class T3>                        \
  A(const T1 &x1, const T2 &x2, const T3 &x3) : B(x1, x2, x3) {} \
  template <class T1, class T2, class T3, class T4>              \
  A(const T1 &x1, const T2 &x2, const T3 &x3, const T4 &x4)      \
      : B(x1, x2, x3, x4) {}

#endif  // HAVE_CONFIG_HPP
