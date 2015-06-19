/** \file
  \brief  Convert vector/matrix-Types to double SEXP types 
*/

double asDouble(int x){return double(x);}
double asDouble(double x){return x;}
double asDouble(AD<double> x){return CppAD::Value(x);}
double asDouble(AD<AD<double> > x){return CppAD::Value(CppAD::Value(x));}
double asDouble(AD<AD<AD<double> > > x){return CppAD::Value(CppAD::Value(CppAD::Value(x)));}

/** \brief Convert TMB matrix, vector, scalar or int to R style */
template<class Type>
SEXP asSEXP(const matrix<Type> &a) 
{
   int nr=a.rows();
   int nc=a.cols();
   int size = nr * nc;
   SEXP val;
   PROTECT(val = allocVector(REALSXP,size));
   double *p = REAL(val);
   for(int i=0;i<nr;i++)
     for(int j=0;j<nc;j++)
       p[i+j*nr]=asDouble(a(i,j));
   SEXP dim;
   PROTECT(dim = allocVector(INTSXP,2));
   INTEGER(dim)[0] = nr; INTEGER(dim)[1] = nc;
   setAttrib(val, R_DimSymbol, dim);
   UNPROTECT(2);
   return val;
}

// Report vector of numeric types: Make R-vector
#define asSEXP_VECTOR_OF_NUMERIC(Type)			\
SEXP asSEXP(const vector<Type> &a)			\
{							\
  int size = a.size();					\
  SEXP val;						\
  PROTECT(val = allocVector(REALSXP,size));		\
  double *p = REAL(val);				\
  for (int i = 0; i < size; i++) p[i] = asDouble(a[i]);	\
  UNPROTECT(1);						\
  return val;						\
}
asSEXP_VECTOR_OF_NUMERIC(int)
asSEXP_VECTOR_OF_NUMERIC(double)
template<class Type>
asSEXP_VECTOR_OF_NUMERIC(AD<Type>)
#undef asSEXP_VECTOR_OF_NUMERIC
// Report vector of anything else: Make R-list
template<class Type>
SEXP asSEXP(const vector<Type> &a)
{
   int size = a.size();
   SEXP val;
   PROTECT(val = allocVector(VECSXP, size));
   for (int i = 0; i < size; i++) SET_VECTOR_ELT(val, i, asSEXP(a[i]));
   UNPROTECT(1);
   return val;
}

SEXP asSEXP(const double &a)
{
   SEXP val;
   PROTECT(val=allocVector(REALSXP,1));
   REAL(val)[0]=a;
   UNPROTECT(1);
   return val;
}
SEXP asSEXP(const int &a)
{
   SEXP val;
   PROTECT(val=allocVector(INTSXP,1));
   INTEGER(val)[0]=a;
   UNPROTECT(1);
   return val;
}
// EXPERIMENT
template<class Type>
SEXP asSEXP(const AD<Type> &a){
  return asSEXP(CppAD::Value(a));
}

/** \brief Construct c++-vector from SEXP object */
template <class Type>
vector<Type> asVector(SEXP x)
{
   if(!isReal(x))error("NOT A VECTOR!");
   int n=length(x);
   vector<Type> y(n);
   for(int i=0;i<n;i++) y[i]=Type(REAL(x)[i]);
   return y;
}

/** \brief Vector <-> Matrix conversion (for row-major matrices) */
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
/** \brief Construct c++-matrix from SEXP object */
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
   PROTECT(val = allocVector(REALSXP,a.size()));
   double *p = REAL(val);
   for(int i=0;i<a.size();i++)p[i]=asDouble(a[i]);
   SEXP dim;
   PROTECT(dim = allocVector(INTSXP,a.dim.size()));
   for(int i=0;i<a.dim.size();i++)INTEGER(dim)[i]=a.dim[i];
   setAttrib(val, R_DimSymbol, dim);
   UNPROTECT(2);
   return val;
}

/** Create R-triplet sparse matrix from Eigen sparse matrix */
template<class Type>
SEXP asSEXP(Eigen::SparseMatrix<Type> x){
  typedef typename Eigen::SparseMatrix<Type>::InnerIterator Iterator;
  // Allocate return object
  int nnz = x.nonZeros();
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
  int k=0;
  for (int cx=0; cx<x.outerSize(); cx++)
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
