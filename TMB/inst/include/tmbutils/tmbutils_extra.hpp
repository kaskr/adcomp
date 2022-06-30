// Copyright (C) 2013-2015 Kasper Kristensen
// License: GPL-2

namespace tmbutils {

template<class Type, class T1, class T2>
vector<Type> dnorm(vector<Type> x, T1 mean, T2 sd, int give_log=0)
{
  vector<Type> logres;
  x=(x-mean)/sd;
  logres=-log(Type(sqrt(2*M_PI))*sd)-Type(.5)*x*x;
  if(give_log)return logres; else return exp(logres);
}

} // End namespace

// Convenience utilites
#include "spmat.hpp"
#include "kronecker.hpp"
#include "matexp.hpp"
#include "order.hpp"
#include "splines.hpp"
#include "interpol.hpp"
#include "density.hpp"
#include "romberg.hpp"
#include "autodiff.hpp"
