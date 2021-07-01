#ifndef TINY_AD_TWEEDIE_H
#define TINY_AD_TWEEDIE_H

#include "../gamma/gamma.hpp"

namespace tweedie_utils {
using gamma_utils::lgammafn;

template<class S, class T>
T fmax2(S x, T y) {return (x < y) ? y : x;}
template<class S>
int imax2(S x, int y) {return (x < y) ? y : x;}
template<class S>
int imin2(S x, int y) {return (x < y) ? x : y;}

#include "tweedie.cpp"
#undef TWEEDIE_DROP
#undef TWEEDIE_INCRE
#undef TWEEDIE_NTERM

} // End namespace tweedie_utils

#endif
