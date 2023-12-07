/** \file
    \brief Include this file to extract declarations, definitions and selected code for pre-compilation
*/
#undef WITH_LIBTMB
#undef TMB_PRECOMPILE_ATOMICS
#undef CSKIP_ATOMIC
#undef IF_TMB_PRECOMPILE_ATOMICS
#undef TMB_EXTERN
// Redefine
#undef  WITH_LIBTMB
#define TMB_PRECOMPILE_ATOMICS
#undef  HAVE_PRECOMPILED_ATOMICS
#define CSKIP_ATOMIC(...) __VA_ARGS__
#define IF_TMB_PRECOMPILE_ATOMICS(...) __VA_ARGS__
#define TMB_EXTERN
