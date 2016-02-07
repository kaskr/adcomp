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
#define CSKIP(x) ;
#define TMB_EXTERN extern
#else
#define CSKIP(x) x
#define TMB_EXTERN
#endif

/* Early inclusion of Rprintf and REprintf */
#include <R_ext/Print.h>
#include "Rstream.hpp"

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
				  abort();}
#else
#undef NDEBUG
#define NDEBUG 1
#endif
#include <Eigen/Dense>
#include <Eigen/Sparse>

/* Include the CppAD library. (Always turn off debug for cppad) */
#undef NDEBUG
#define NDEBUG 1
#include "cppad/cppad.hpp"

/* Include the R library _after_ Eigen and CppAD. Otherwise, the R
   macros can cause conflicts (as they do not respect the Eigen and
   CppAD namespace limits). E.g., the 'length' macro conflicts with
   CppAD when compiling with '-std=c++11'. */
#include <R.h>
#include <Rinternals.h>
void eigen_REprintf(const char* x)CSKIP({REprintf(x);})

#include "tmbutils/tmbutils.hpp"
using tmbutils::matrix;
using tmbutils::vector;
using CppAD::AD;
using CppAD::ADFun;
namespace CppAD{
  /* Add to CppAD so that 'Variable' works for any 'Type' */
  bool Variable(double x)CSKIP({ return false; })
}
#include "convert.hpp" // asSEXP, asMatrix, asVector
#include "config.hpp"
#include "atomic_math.hpp"
#include "expm.hpp"
#include "Vectorize.hpp"
#include "dnorm.hpp"   // harmless
#include "lgamma.hpp"  // harmless
#include "start_parallel.hpp"
#include "tmb_core.hpp"
#include "convenience.hpp"
#include "distributions_R.hpp"
#include "tmbutils/tmbutils_extra.hpp"
#include "tmbutils/R_inla.hpp"
#include "precompile.hpp" // Must come last
using tmbutils::array;
using Eigen::Matrix;
using Eigen::Array;

/* Cleanup  */
#undef CSKIP        // Nothing more to precompile. Must disable
#define CSKIP(x) x  // to not confuse REGISTER_ATOMIC etc.
