// using boost::numeric::ublas::vector;
// using boost::numeric::ublas::matrix;

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

