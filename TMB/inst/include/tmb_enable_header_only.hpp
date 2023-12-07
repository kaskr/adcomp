/** \file
    \brief Include this file to extract declarations only
*/
#undef WITH_LIBTMB
#undef TMB_PRECOMPILE_ATOMICS
#undef CSKIP_ATOMIC
#undef IF_TMB_PRECOMPILE_ATOMICS
#undef TMB_EXTERN
// Redefine
#define WITH_LIBTMB
#undef  TMB_PRECOMPILE_ATOMICS
#define HAVE_PRECOMPILED_ATOMICS
#define CSKIP_ATOMIC(...) ;
#define IF_TMB_PRECOMPILE_ATOMICS(...)
#define TMB_EXTERN extern
