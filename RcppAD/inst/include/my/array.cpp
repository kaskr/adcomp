/* Extend existing vector class with dim member and modify operator() method */
template<class Type>
struct array:Map< Array<Type,Dynamic,1> >{

  typedef Array<Type,Dynamic,1> Base;
  typedef Map< Base > MapBase;

  vector<int> dim;
  vector<int> mult;

  Base vectorcopy; /* Used by the copy-constructor to keep a protected copy */

  //array(Base x, vector<int> dim_):Base::vector(x){dim=dim_;}
  void setdim(vector<int> dim_){
    dim=dim_;
    mult.resize(dim.size());
    mult[dim.size()-1]=1;
    for(int k=dim.size()-1;k>0;k--){
      mult[k-1]=mult[k]*dim[k];
    }
  }
  void initZeroArray(vector<int> dim_){
    vectorcopy.resize(dim_.prod());
    vectorcopy.setZero();
    new (this) MapBase(&vectorcopy[0],vectorcopy.size()); /* Eigen manual: Despite appearances, this does not invoke the memory allocator... */
    setdim(dim_);
  }
  /* Default constructor: e.g. array<double> a; */
  array(){};

  /* Constructor from dimension alone */
  array(vector<int> dim_):MapBase(NULL,0){
    initZeroArray(dim_);
  }
  array(int n1):MapBase(NULL,0){
    vector<int> dim_(1);
    dim_ << n1;
    initZeroArray(dim_);
  }
  array(int n1, int n2):MapBase(NULL,0){
    vector<int> dim_(2);
    dim_ << n1,n2;
    initZeroArray(dim_);
  }
  array(int n1, int n2, int n3):MapBase(NULL,0){
    vector<int> dim_(3);
    dim_ << n1,n2,n3;
    initZeroArray(dim_);
  }
  array(int n1, int n2, int n3, int n4):MapBase(NULL,0){
    vector<int> dim_(4);
    dim_ << n1,n2,n3,n4;
    initZeroArray(dim_);
  }

  /* Default construction: always deep copy ! */
  template<class T>
  array(T &x, vector<int> dim_):vectorcopy(x),MapBase(NULL,0){
    if(x.size()>0){
      new (this) MapBase(&vectorcopy[0],x.size()); /* Eigen manual: Despite appearances, this does not invoke the memory allocator... */
    }
    setdim(dim_);
  }
  /* Deep copy existing array - as above with dim_=x.dim  */
  template<class T>
  array(array<T> &x):vectorcopy(x),MapBase(NULL,0){
    if(x.size()>0){
      new (this) MapBase(&vectorcopy[0],x.size()); /* Eigen manual: Despite appearances, this does not invoke the memory allocator... */
    }
    setdim(x.dim);
  }
  /* Sometimes we want a reference only... See operator() */
  array(Type *p, vector<int> dim_):MapBase(p,dim_.prod()){
    setdim(dim_);
  }
  
  void print(){
    std::cout << "Array dim: ";
    for(int i=0;i<dim.size();i++)std::cout << dim[i] << " ";
    std::cout << "\n";
    std::cout << "Array val: ";
    for(int i=0;i<this->MapBase::size();i++)std::cout << this->MapBase::operator[](i) << " ";
    std::cout << "\n";
  };
  int rows(){
    return dim[0];
  }

  MapBase row(int i){
    int nslice=this->MapBase::size()/dim[0];
    return MapBase(&(this->MapBase::operator()(i*nslice)),nslice);
  }
  array<Type> operator()(int i){
    MapBase v=this->row(i);
    vector<int> newdim;
    if(dim.size()>1){
      newdim=dim.segment(1,dim.size()-1);
    } else {
      newdim.resize(1);
      newdim << 1;
    }
    return array(&v[0],newdim);
  }
  array<Type> operator()(int i1, int i2){
    return (*this)(i1)(i2);
  }
  array<Type> operator()(int i1, int i2, int i3){
    return (*this)(i1)(i2)(i3);
  }
  array<Type> operator()(int i1, int i2, int i3, int i4){
    return (*this)(i1)(i2)(i3)(i4);
  }
  int index(vector<int> tup){
    return (tup*mult).sum();
  }
  vector<int> tuple(int i){
    vector<int> tup(dim.size());
    for(int k=0;k<dim.size();k++){
      tup[k]=i/mult[k];
      i=i-tup[k]*mult[k];
    }    
    return tup;
  }

  array<Type> perm(vector<int> p){
    vector<Type> x(this->size());
    array<Type> ans(x,dim(p));       /* Create new array with permuted dimension */
    for(int i=0;i<this->size();i++){ /* Loop through values of old array */
      ans[ans.index(tuple(i)(p))]=this->operator[](i);
    }
    return ans;
  }
  /* Special case: array transpose */
  array<Type> transpose(){
    vector<int> p(dim.size());
    for(int i=0;i<p.size();i++)p[i]=i;
    return this->perm(p.reverse());
  }
  /* Other special case: array rotate */
  int mod(int i,int n){return ((i%n)+n)%n;}
  array<Type> rotate(int n){
    vector<int> p(dim.size());
    for(int i=0;i<p.size();i++)p[i]=mod(i-n,p.size());
    return this->perm(p);
  }

#define INHERIT(OP)					\
  template <class T>					\
  array<Type> OP(T y){return array(MapBase::OP(y),dim);}
  INHERIT(operator+)
  INHERIT(operator-)
  INHERIT(operator*)
  INHERIT(operator/)
  INHERIT(operator=)
#undef INHERIT

};

