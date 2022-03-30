/** \file
    \brief Include this file to extract declarations, definitions and selected code for pre-compilation
*/
#undef WITH_LIBTMB
#undef TMB_PRECOMPILE
#undef CSKIP
#undef IF_TMB_PRECOMPILE
#undef TMB_EXTERN
// Redefine
#undef  WITH_LIBTMB
#define TMB_PRECOMPILE
#define CSKIP(...) __VA_ARGS__
#define IF_TMB_PRECOMPILE(...) __VA_ARGS__
#define TMB_EXTERN
