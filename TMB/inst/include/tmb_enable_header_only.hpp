/** \file
    \brief Include this file to extract declarations only
*/
#undef WITH_LIBTMB
#undef TMB_PRECOMPILE
#undef CSKIP
#undef IF_TMB_PRECOMPILE
#undef TMB_EXTERN
// Redefine
#define WITH_LIBTMB
#undef  TMB_PRECOMPILE
#define CSKIP(...) ;
#define IF_TMB_PRECOMPILE(...)
#define TMB_EXTERN extern
