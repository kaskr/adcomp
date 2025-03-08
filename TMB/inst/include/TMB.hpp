// Copyright (C) 2013-2015 Kasper Kristensen
// License: GPL-2
/* Utility: Compile time test for Type=double */
template<class Type>
struct isDouble{
  enum{value=false};
};
template<>
struct isDouble<double>{
  enum{value=true};
};

/** \file 
* \brief Includes and sets all stuff needed to compile the user defined objective function.
*/
/* To be removed */
#define TMB_DEBUG 0
#define TMB_PRINT(x)std::cout << #x << ": " << x << "\n"; std::cout.flush();

/* Conditionally skip compilation */
#ifdef WITH_LIBTMB
#define CSKIP(...) ;
#define TMB_EXTERN extern
#else
#define CSKIP(...) __VA_ARGS__
#define TMB_EXTERN
#endif
#ifdef TMB_PRECOMPILE_ATOMICS
#define IF_TMB_PRECOMPILE_ATOMICS(...) __VA_ARGS__
#else
#define IF_TMB_PRECOMPILE_ATOMICS(...)
#endif
#ifdef HAVE_PRECOMPILED_ATOMICS
#define CSKIP_ATOMIC(...) ;
#else
#define CSKIP_ATOMIC(...) __VA_ARGS__
#endif

/* Must come before Rinternals.h */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Early inclusion of Rprintf and REprintf */
#include <R_ext/Print.h>
#include "Rstream.hpp"

/* Flag to bypass abort() */
#ifndef TMB_ABORT
#define TMB_ABORT abort()
#endif

/* Include the Eigen library. */
#ifdef TMB_SAFEBOUNDS
#undef NDEBUG
#undef eigen_assert
void eigen_REprintf(const char* x);
#define eigen_assert(x) if (!(x)) { eigen_REprintf("TMB has received an error from Eigen. "); \
                                  eigen_REprintf("The following condition was not met:\n");          \
                                  eigen_REprintf(#x);                                                \
                                  eigen_REprintf("\nPlease check your matrix-vector bounds etc., "); \
                                  eigen_REprintf("or run your program through a debugger.\n");       \
				  TMB_ABORT;}
#define TMBAD_ASSERT2(x,msg)                                            \
if (!(x)) {                                                             \
  Rcerr << "TMBad assertion failed.\n";                                 \
  Rcerr << "The following condition was not met: " << #x << "\n";       \
  Rcerr << "Possible reason: " msg << "\n";                             \
  Rcerr << "For more info run your program through a debugger.\n";      \
  TMB_ABORT;                                                            \
}
#define TMBAD_ASSERT(x) TMBAD_ASSERT2(x,"Unknown")
#else
#undef NDEBUG
#define NDEBUG 1
#define TMBAD_ASSERT2(x,msg) (void) (x);
#define TMBAD_ASSERT(x) (void) (x);
#endif
/* Provide access to file 'DisableStupidWarnings.h' which has been
   patched by RcppEigen to satisfy CRAN policy. This file may need
   regular updating. The renaming is to aviod a CRAN note. */
#ifdef TMB_EIGEN_DISABLE_WARNINGS
#ifndef EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS 1
#endif
#include "EigenWarnings/DisableStupidWarnings"
#endif
/* We cannot use Eigen's parallel matrix multiply for AD types (GH390). */
#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif
#include <Eigen/Dense>

// Default: Include Eigen/Sparse normally
#ifndef TMB_SPARSE_STORAGE_INDEX
#include <Eigen/Sparse>
#else
// Alternative: Include Eigen/Sparse with custom sparse matrix integer type
#define SparseMatrix SparseMatrix_rename
#include <Eigen/Sparse>
#undef SparseMatrix
namespace Eigen {
template<class T, int Flags = 0, class StorageIndex = TMB_SPARSE_STORAGE_INDEX>
using SparseMatrix = SparseMatrix_rename<T, Flags, StorageIndex>;
}
#endif

/* Workaround side effect when -DEIGEN_USE_LAPACKE is set */
#undef I

/* Select AD framework: TMBAD or CPPAD  */
#ifndef CPPAD_FRAMEWORK
#ifndef TMBAD_FRAMEWORK
#define CPPAD_FRAMEWORK
#endif
#endif

// TMBad config variables
/** \brief Use deterministic hash codes for tape optimizer ? */
TMB_EXTERN bool tmbad_deterministic_hash;

/* Include the CppAD library. (Always turn off debug for cppad) */
#undef NDEBUG
#define NDEBUG 1
#include "cppad/cppad.hpp"
#ifdef TMBAD_FRAMEWORK
#define TMBAD_DETERMINISTIC_HASH tmbad_deterministic_hash
#include "TMBad/TMBad.hpp"
#include "TMBad/tmbad_allow_comparison.hpp"
#include "TMBad/eigen_numtraits.hpp"
#undef error
#include "TMBad/compile.hpp"
#include "TMBad/graph2dot.hpp"
#include "TMBad/compression.hpp"
#include "TMBad/ad_blas.hpp"
#ifndef WITH_LIBTMB
#include "TMBad/TMBad.cpp"
#endif
#define error Rf_error
// Workaround to make CppAD::Integer working with TMBad
namespace CppAD {
int Integer(const TMBad::ad_aug &x) CSKIP ({ return (int) x.Value(); })
TMBad::ad_aug abs(const TMBad::ad_aug &x) CSKIP ({ return TMBad::fabs(x); })
#define TMBAD_CONDEXP(NAME)                             \
TMBad::ad_aug CondExp ## NAME(                          \
  const TMBad::ad_aug &x0,                              \
  const TMBad::ad_aug &x1,                              \
  const TMBad::ad_aug &x2,                              \
  const TMBad::ad_aug &x3) CSKIP ( {                    \
      return TMBad::CondExp ## NAME(x0, x1, x2, x3);    \
})
TMBAD_CONDEXP(Eq)
TMBAD_CONDEXP(Ne)
TMBAD_CONDEXP(Lt)
TMBAD_CONDEXP(Gt)
TMBAD_CONDEXP(Le)
TMBAD_CONDEXP(Ge)
#undef TMBAD_CONDEXP
bool Variable(const TMBad::ad_aug &x) CSKIP ({ return !x.constant(); })
}
// FIXME: Move to TMBad source?
namespace TMBad {
  /* Add 'isfinite', 'isinf' and 'isnan' to TMBad */
  using std::isfinite;
  bool isfinite(const TMBad::ad_aug &x)CSKIP({ return isfinite(x.Value()); })
  using std::isinf;
  bool isinf(const TMBad::ad_aug &x)CSKIP({ return isinf(x.Value()); })
  using std::isnan;
  bool isnan(const TMBad::ad_aug &x)CSKIP({ return isnan(x.Value()); })
}
// Add missing numeric_limits specialization
namespace std {
template<>
struct numeric_limits<TMBad::ad_aug> : numeric_limits<TMBad::Scalar> { };
}
#endif

/* Include the R library _after_ Eigen and CppAD. Otherwise, the R
   macros can cause conflicts (as they do not respect the Eigen and
   CppAD namespace limits). E.g., the 'length' macro conflicts with
   CppAD when compiling with '-std=c++11'. */
#include <R.h>
#include <Rinternals.h>
/* See WRE manual */
#include <Rversion.h>
#if R_VERSION < R_Version(4, 5, 0)
#define R_ParentEnv(x) ENCLOS(x)
#endif

#include "toggle_thread_safe_R.hpp"
void eigen_REprintf(const char* x)CSKIP({REprintf("%s",x);})

#include "tmbutils/tmbutils.hpp"
#include "tmbutils/vectorize.hpp"
using tmbutils::matrix;
using tmbutils::vector;
using CppAD::AD;
using CppAD::ADFun;
namespace CppAD{
  /* Add to CppAD so that 'Variable' works for any 'Type' */
  bool Variable(double x)CSKIP({ return false; })
  /* Add 'isfinite', 'isinf' and 'isnan' to CppAD */
  using std::isfinite;
  template <class T>
  bool isfinite(const AD<T> &x)CSKIP({ return isfinite(Value(x)); })
  using std::isinf;
  template <class T>
  bool isinf(const AD<T> &x)CSKIP({ return isinf(Value(x)); })
  using std::isnan;
  template <class T>
  bool isnan(const AD<T> &x)CSKIP({ return isnan(Value(x)); })
}
#include "convert.hpp" // asSEXP, asMatrix, asVector
#include "config.hpp"
#include "tmbutils/getListElement.hpp"
#include "atomic_math.hpp"
#include "expm.hpp"
#include "atomic_convolve.hpp"
#include "tiny_ad/atomic.hpp"
#include "tiny_ad/integrate/integrate.hpp"
#include "dynamic_data.hpp" // Requires atomic namespace
#include "Vectorize.hpp"
#include "dnorm.hpp"   // harmless
#include "lgamma.hpp"  // harmless
#include "start_parallel.hpp"
#include "tmbutils/newton.hpp" // Newton solver + Laplace used by TransformADFunObject
#ifndef TMB_SKINNY
#include "tmb_core.hpp"
#endif
#include "distributions_R.hpp"
#include "convenience.hpp"    // Requires besselK
#include "tmbutils/tmbutils_extra.hpp"
#include "tmbutils/R_inla.hpp"
#include "tmbutils/sparse_matrix_exponential.hpp"
#include "tmbutils/concat.hpp"
#include "precompile.hpp" // Must come last
using tmbutils::array;
using Eigen::Matrix;
using Eigen::Array;

/* Cleanup  */

// Nothing more to precompile
#undef CSKIP
#define CSKIP(...) __VA_ARGS__
#undef CSKIP_ATOMIC
#define CSKIP_ATOMIC(...) __VA_ARGS__
