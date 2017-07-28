// Copyright (C) 2017 Bruce Swihart
// License: GPL-2

/** \file
    \brief Univariate bridge density (Wang & Louis (2003))
    \ingroup R_style_distribution
*/
template<class Type>
Type dbridge(Type x, Type phi, int give_log=0)
{
  Type logres;
  logres=-log(Type(2*M_PI))+log(sin(Type(phi*M_PI)))-log(cosh(Type(phi*x))+cos(Type(phi*M_PI)));
  //return 1/(2*M_PI)*sin(phi*M_PI)/(cosh(phi*x)+cos(phi*M_PI));
  if(give_log)return logres; else return exp(logres);
}
VECTORIZE3_tti(dbridge)
