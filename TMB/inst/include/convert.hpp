// Copyright (C) 2013-2015 Kasper Kristensen
// License: GPL-2

/** \file
  \brief  Convert vector/matrix-Types to double SEXP types 
*/

#ifdef WITH_LIBTMB
double asDouble(int x);
double asDouble(double x);
double asDouble(AD<double> x);
double asDouble(AD<AD<double> > x);
double asDouble(AD<AD<AD<double> > > x);
#else
double asDouble(int x){return double(x);}
double asDouble(double x){return x;}
double asDouble(AD<double> x){return CppAD::Value(x);}
double asDouble(AD<AD<double> > x){return CppAD::Value(CppAD::Value(x));}
double asDouble(AD<AD<AD<double> > > x){return CppAD::Value(CppAD::Value(CppAD::Value(x)));}
#endif

/** \brief Convert TMB matrix, vector, scalar or int to R style */
template<class Type>
SEXP asSEXP(const matrix<Type> &a) 
{
   R_xlen_t nr = a.rows();
   R_xlen_t nc = a.cols();
   SEXP val;
   PROTECT(val = allocMatrix(REALSXP, nr, nc));
   double *p = REAL(val);
   for(R_xlen_t i=0; i<nr; i++)
     for(R_xlen_t j=0; j<nc; j++)
       p[i + j * nr] = asDouble(a(i,j));
   UNPROTECT(1);
   return val;
}

// Report vector of numeric types: Make R-vector
#define asSEXP_VECTOR_OF_NUMERIC(Type)          \
SEXP asSEXP(const vector<Type> &a) CSKIP(       \
{                                               \
  R_xlen_t size = a.size();                     \
  SEXP val;                                     \
  PROTECT(val = allocVector(REALSXP,size));     \
  double *p = REAL(val);                        \
  for (R_xlen_t i = 0; i < size; i++)           \
    p[i] = asDouble(a[i]);                      \
  UNPROTECT(1);                                 \
  return val;                                   \
})
asSEXP_VECTOR_OF_NUMERIC(int)
asSEXP_VECTOR_OF_NUMERIC(double)
template<class Type>
asSEXP_VECTOR_OF_NUMERIC(AD<Type>)
#undef asSEXP_VECTOR_OF_NUMERIC
// Report vector of anything else: Make R-list
template<class Type>
SEXP asSEXP(const vector<Type> &a)
{
   R_xlen_t size = a.size();
   SEXP val;
   PROTECT(val = allocVector(VECSXP, size));
   for (R_xlen_t i = 0; i < size; i++)
     SET_VECTOR_ELT(val, i, asSEXP(a[i]));
   UNPROTECT(1);
   return val;
}

SEXP asSEXP(const double &a) CSKIP(
{
   SEXP val;
   PROTECT(val=allocVector(REALSXP,1));
   REAL(val)[0]=a;
   UNPROTECT(1);
   return val;
})
SEXP asSEXP(const int &a) CSKIP(
{
   SEXP val;
   PROTECT(val=allocVector(INTSXP,1));
   INTEGER(val)[0]=a;
   UNPROTECT(1);
   return val;
})
// EXPERIMENT
template<class Type>
SEXP asSEXP(const AD<Type> &a){
  return asSEXP(CppAD::Value(a));
}

/** \brief Construct c++-vector from SEXP object */
template <class Type>
vector<Type> asVector(SEXP x)
{
   if(!isReal(x)) error("NOT A VECTOR!");
   R_xlen_t n = XLENGTH(x);
   typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > MapVector;
   MapVector tmp(REAL(x), n);
   vector<Type> y = tmp.cast<Type>();
   return y;
}

/** \brief Vector <-> Matrix conversion (for row-major matrices) */
template<class Type>
matrix<Type> asMatrix(const vector<Type> &x, int nr, int nc)
{
  matrix<Type> xm = x.matrix();
  xm.resize(nr, nc);
  return xm;
}

/** \brief Construct C++-matrix from SEXP object */
template <class Type>
matrix<Type> asMatrix(SEXP x)
{
   if (!isMatrix(x))
     error("x must be a matrix in 'asMatrix(x)'");
   R_xlen_t nr = nrows(x); // nrows is int
   R_xlen_t nc = ncols(x); // ncols is int
   matrix<Type> y(nr, nc);
   for(R_xlen_t i=0; i<nr; i++)
     for(R_xlen_t j=0; j<nc; j++)
       y(i, j) = Type(REAL(x)[i + nr * j]);
   return y;
}

template<class Type>
SEXP asSEXP(const tmbutils::array<Type> &a)
{
   SEXP val;
   PROTECT(val = allocVector(REALSXP, a.size()));
   double *p = REAL(val);
   for(R_xlen_t i=0; i<a.size(); i++)
     p[i] = asDouble(a[i]);
   SEXP dim;
   PROTECT(dim = allocVector(INTSXP, a.dim.size()));
   for(int i=0; i<a.dim.size(); i++)
     INTEGER(dim)[i] = a.dim[i];
   setAttrib(val, R_DimSymbol, dim);
   UNPROTECT(2);
   return val;
}

/** Create R-triplet sparse matrix from Eigen sparse matrix */
template<class Type>
SEXP asSEXP(Eigen::SparseMatrix<Type> x){
  typedef typename Eigen::SparseMatrix<Type>::InnerIterator Iterator;
  // Allocate return object
  R_xlen_t nnz = x.nonZeros();
  SEXP ans = PROTECT(R_do_new_object(R_do_MAKE_CLASS("dgTMatrix")));
  SEXP dim = PROTECT(allocVector(INTSXP, 2));
  SEXP dimnames = PROTECT(allocVector(VECSXP, 2));
  SEXP values = PROTECT(allocVector(REALSXP, nnz));
  SEXP i = PROTECT(allocVector(INTSXP, nnz));
  SEXP j = PROTECT(allocVector(INTSXP, nnz));
  SEXP factors = PROTECT(allocVector(VECSXP, 0));
  R_do_slot_assign(ans, install("i"), i);
  R_do_slot_assign(ans, install("j"), j);
  R_do_slot_assign(ans, install("Dim"), dim);
  R_do_slot_assign(ans, install("Dimnames"), dimnames);
  R_do_slot_assign(ans, install("x"), values);
  R_do_slot_assign(ans, install("factors"), factors);
  // Insert
  INTEGER(dim)[0] = x.rows();
  INTEGER(dim)[1] = x.cols();
  R_xlen_t k = 0;
  for (R_xlen_t cx=0; cx<x.outerSize(); cx++)
    {
      for (Iterator itx(x,cx); itx; ++itx)
	{
	  INTEGER(i)[k] = itx.row();
	  INTEGER(j)[k] = itx.col();
	  REAL(values)[k] = asDouble(itx.value());
	  k++;
	}
    }
  UNPROTECT(7);
  return ans;
}


/** \brief Convert C++ object to R object and protect the result from garbage collection
    \note Each call to this function must be followed by UNPROTECT(1)
*/
template<class T>
SEXP asSEXP_protect(const T &x) {
  SEXP ans;
  PROTECT( ans = asSEXP(x) );
  return ans;
}
