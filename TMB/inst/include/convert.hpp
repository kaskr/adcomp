/* Convert vector/matrix-Types to double SEXP types */
double asDouble(int x){return double(x);}
double asDouble(double x){return x;}
double asDouble(AD<double> x){return CppAD::Value(x);}
double asDouble(AD<AD<double> > x){return CppAD::Value(CppAD::Value(x));}
double asDouble(AD<AD<AD<double> > > x){return CppAD::Value(CppAD::Value(CppAD::Value(x)));}
template<class Type>
SEXP asSEXP(const matrix<Type> &a) 
{
   int nr=a.rows();
   int nc=a.cols();
   int size = nr * nc;
   SEXP val;
   PROTECT(val = NEW_NUMERIC(size));
   double *p = NUMERIC_POINTER(val);
   for(int i=0;i<nr;i++)
     for(int j=0;j<nc;j++)
       p[i+j*nr]=asDouble(a(i,j));
   SEXP dim;
   PROTECT(dim = NEW_INTEGER(2));
   INTEGER(dim)[0] = nr; INTEGER(dim)[1] = nc;
   SET_DIM(val, dim);
   UNPROTECT(2);
   return val;
}
template<class Type>
SEXP asSEXP(const vector<Type> &a) 
{
  int size = a.size();
  SEXP val;
  PROTECT(val = NEW_NUMERIC(size));
  double *p = NUMERIC_POINTER(val);
  for (int i = 0; i < size; i++) p[i] = asDouble(a[i]);
  SEXP len;
  PROTECT(len = NEW_INTEGER(1));
  INTEGER(len)[0] = size;
  SET_LENGTH(val, size);
  UNPROTECT(2); 
  return val;
}
template<class Type>
SEXP asSEXP(const Type &a)
{
   SEXP val;
   PROTECT(val=NEW_NUMERIC(1));
   REAL(val)[0]=asDouble(a);
   UNPROTECT(1);
   return val;
}
SEXP asSEXP(const int &a)
{
   SEXP val;
   PROTECT(val=NEW_INTEGER(1));
   INTEGER(val)[0]=a;
   UNPROTECT(1);
   return val;
}
// EXPERIMENT
template<class Type>
SEXP asSEXP(const AD<Type> &a){
  return asSEXP(CppAD::Value(a));
}
template<template<class> class Vector, class Type>
SEXP asSEXP(const Vector<Type> &a)
{
   int size = a.size();
   SEXP val;
   PROTECT(val = NEW_NUMERIC(size));
   double *p = NUMERIC_POINTER(val);
   for (int i = 0; i < size; i++) p[i] = asDouble(a[i]);
   SEXP len;
   PROTECT(len = NEW_INTEGER(1));
   INTEGER(len)[0] = size;
   SET_LENGTH(val, size);
   UNPROTECT(2); 
   return val;
}

/* Construct c++-vector from SEXP object */
template <class Type>
vector<Type> asVector(SEXP x)
{
   if(!isReal(x))error("NOT A VECTOR!");
   int n=length(x);
   vector<Type> y(n);
   for(int i=0;i<n;i++) y[i]=Type(REAL(x)[i]);
   return y;
}

/* Vector <-> Matrix conversion (for row-major matrices) */
template<class Type>
matrix<Type> asMatrix(const vector<Type> &x, int nr, int nc)
{
  if(nr*nc!=x.size())error("nr*nc!=n in asMatrix");
  matrix<Type> res(nr,nc);
  for(int i=0;i<nr;i++)
    for(int j=0;j<nc;j++)
      res(i,j)=x[i*nc+j];
  return res;
}

// kasper: MOVE TO asMatrix.hpp ?
/* Construct c++-matrix from SEXP object */
template <class Type>
matrix<Type> asMatrix(SEXP x)
{
   if(!isMatrix(x))error("NOT A MATRIX!");
   int nr=nrows(x);
   int nc=ncols(x);
   matrix<Type> y(nr,nc);
   for(int i=0;i<nr;i++)
     for(int j=0;j<nc;j++)
       y(i,j)=Type(REAL(x)[i+nr*j]);
   return y;
}


template<class Type>
vector<Type> asVector(matrix<Type> x)
{
  int nr=x.size1();
  int nc=x.size2();
  vector<Type> res(nr*nc);
  for(int i=0;i<nr;i++)
    for(int j=0;j<nc;j++)
      res[i*nc+j]=x(i,j);
  return res;
}

template<class Type>
SEXP asSEXP(const tmbutils::array<Type> &a)
{
   SEXP val;
   PROTECT(val = NEW_NUMERIC(a.size()));
   double *p = REAL(val);
   for(int i=0;i<a.size();i++)p[i]=asDouble(a[i]);
   SEXP dim;
   PROTECT(dim = NEW_INTEGER(a.dim.size()));
   for(int i=0;i<a.dim.size();i++)INTEGER(dim)[i]=a.dim[i];
   SET_DIM(val, dim);
   UNPROTECT(2);
   return val;
}
