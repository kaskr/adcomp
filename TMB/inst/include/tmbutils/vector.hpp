// Copyright (C) 2013-2015 Kasper Kristensen
// License: GPL-2

/** \file 
   \brief Defines TMB vectors 
*/

using namespace Eigen;

/** \brief Vector class used by TMB.

    The TMB vector class is implemented as an Eigen Array of
    dynamic length. In particular, vectorized operations are inherited
    from the Eigen library.
*/
template <class Type>
struct vector : Array<Type,Dynamic,1>
{
  typedef Type value_type;
  typedef Array<Type,Dynamic,1> Base;
  vector(void):Base() {}

  template<class T1>
  vector(T1 x):Base(x) {}

  template<class T1, class T2>
  vector(T1 x, T2 y):Base(x,y) {}
  
  
  template<class T1>
  vector & operator= (const T1 & other)
  {
    this->Base::operator=(other);
    return *this;
  }

  // /* index-vector subset */
  // template <class T>
  // Type & operator()(T ind){
  //   return this->Base::operator()(ind);
  // }
  using Base::operator();

  // kasper: would be better with references, i.e. allow x(indvec)=y;
  vector<Type> operator()(vector<int> ind){
    vector<Type> ans(ind.size());
    for(int i=0;i<ind.size();i++)ans[i]=this->operator[](ind[i]);
    return ans;
  }

  /* Convert to _other_ vector class
     Examples
     1. vector<Type> to CppAD::vector<Type>
     2. vector<int> to CppAD::vector<double>
   */
  template<template<class> class Vector, class T>
  operator Vector<T>(){
    int n = this->size();
    Vector<T> x(n);
    for(int i=0; i<n; i++) x[i] = T(this->operator[](i));
    return x;
  }
  /* Convert to _this_ vector class
     Examples:
     1. vector<int> to vector<double>
     2. CppAD::vector<Type> to vector<Type> 
     3. CppAD::vector<int> to vector<double> 
  */
  template<template<class> class Vector, class T>
  vector(Vector<T> x):Base(){
    int n = x.size();
    this->resize(n);
    for(int i=0; i<n; i++) this->operator[](i) = Type(x[i]);
  }
};

/** \brief Matrix class used by TMB.

    The TMB matrix class is implemented as an Eigen Matrix of
    dynamic dimension. In particular, linear algebra methods are inherited
    from the Eigen library.
*/
template <class Type>
struct matrix : Matrix<Type,Dynamic,Dynamic>
{
  typedef Matrix<Type,Dynamic,Dynamic> Base;
  matrix(void):Base() {}
  template<class T1>
  matrix(T1 x):Base(x) {}
  template<class T1, class T2>
  matrix(T1 x, T2 y):Base(x,y) {}
  
  template<class T1>
  matrix & operator= (const T1 & other)
  {
    this->Base::operator=(other);
    return *this;
  }

  /**
     The vec operator stacks the matrix columns into a single vector.
  */
  vector<Type> vec(){
    Array<Type,Dynamic,Dynamic> a = this->array();
    a.resize(a.size(), 1);
    return a;
  }
};
