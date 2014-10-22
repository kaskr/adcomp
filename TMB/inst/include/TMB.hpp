/** \file 
* \brief Includes and sets all stuff needed to compile the user defined objective function.
*/
/* To be removed */
#define TMB_DEBUG 0
#define TMB_PRINT(x)std::cout << #x << ": " << x << "\n"; std::cout.flush();

/* Turn on debug for Eigen ? */
#ifdef TMB_SAFEBOUNDS
#undef NDEBUG
#undef eigen_assert
void eigen_Rprintf(const char* x);
#define eigen_assert(x) if (!(x)) { eigen_Rprintf("TMB has received an error from Eigen. "); \
                                  eigen_Rprintf("The following condition was not met:\n");          \
                                  eigen_Rprintf(#x);                                                \
                                  eigen_Rprintf("\nPlease check your matrix-vector bounds etc., "); \
                                  eigen_Rprintf("or run your program through a debugger.\n");       \
				  abort();}
#else
#undef NDEBUG
#define NDEBUG 1
#endif
#include <Eigen/Dense>
#include <Eigen/Sparse>
/* R must come after Eigen because conflict with length macro on mac os x */
#include <R.h>
#include <Rinternals.h>
#include "Rstream.hpp"
void eigen_Rprintf(const char* x){Rprintf(x);}
/* Always turn off debug for cppad */
#undef NDEBUG
#define NDEBUG 1
#include "cppad/cppad.hpp"
#undef NDEBUG
#include "tmbutils/tmbutils.cpp"
using tmbutils::matrix;
using tmbutils::vector;
using CppAD::AD;
using CppAD::ADFun;
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
using tmbutils::array;
using Eigen::Matrix;
using Eigen::Array;
