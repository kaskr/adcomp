#define REGISTER_ATOMIC(NAME)						\
namespace NAME##checkpoint_namespace{					\
  using CppAD::checkpoint;						\
  bool initialized=false;						\
  int m=0;								\
  int n=0;								\
  checkpoint<double> *NAME##_check1;					\
  checkpoint<AD<double> > *NAME##_check2;				\
  checkpoint<AD<AD<double> > > *NAME##_check3;				\
  void first_call(vector<double> x){					\
    vector<double> y = NAME(x);						\
    m=x.size(); n=y.size();						\
  }									\
  template<class Type>							\
  void NAME(const vector<Type> &x, vector<Type> &y){			\
    vector<Type> res=::NAME(x);						\
    for(int i=0;i<res.size();i++)y[i]=res[i];				\
  }									\
  void initialize(vector<double> x){					\
    first_call(x);							\
    vector<AD<double> > x1(m);						\
    vector<AD<double> > y1(n);						\
    for(int i=0;i<m;i++)x1[i]=x[i];					\
    NAME##_check1 = new checkpoint<double>(#NAME "_check1", NAME<AD<double> >, x1, y1);	\
    vector<AD<AD<double> > > x2(m);					\
    vector<AD<AD<double> > > y2(n);					\
    for(int i=0;i<m;i++)x2[i]=x[i];					\
    NAME##_check2 = new checkpoint<AD<double> >(#NAME "_check2", NAME<AD<AD<double> > >, x2, y2); \
    vector<AD<AD<AD<double> > > > x3(m);				\
    vector<AD<AD<AD<double> > > > y3(n);				\
    for(int i=0;i<m;i++)x3[i]=x[i];					\
    NAME##_check3 = new checkpoint<AD<AD<double> > >(#NAME "_check3", NAME<AD<AD<AD<double> > > >, x3, y3); \
  }									\
  vector<double> eval(vector<double> x){				\
    if(!initialized)initialize(x);					\
    vector<double> y(n);						\
    NAME(x,y);								\
    return y;								\
  }									\
  vector<AD<double> > eval(vector<AD<double> > x){			\
    vector<AD<double> > y(n);						\
    NAME##_check1->operator()(x,y);					\
    NAME(x,y);								\
    return y;								\
  }									\
  vector<AD<AD<double> > > eval(vector<AD<AD<double> > > x){		\
    vector<AD<AD<double> > > y(n);					\
    NAME##_check2->operator()(x,y);					\
    NAME(x,y);								\
    return y;								\
  }									\
  vector<AD<AD<AD<double> > > > eval(vector<AD<AD<AD<double> > > > x){	\
    vector<AD<AD<AD<double> > > > y(n);					\
    NAME##_check3->operator()(x,y);					\
    NAME(x,y);								\
    return y;								\
  }									\
};									\
vector<double> NAME(vector<double> x){					\
  return NAME##checkpoint_namespace::eval(x);				\
}									\
vector<AD<double> > NAME(vector<AD<double> > x){			\
  return NAME##checkpoint_namespace::eval(x);				\
}									\
vector<AD<AD<double> > > NAME(vector<AD<AD<double> > > x){		\
  return NAME##checkpoint_namespace::eval(x);				\
}									\
vector<AD<AD<AD<double> > > > NAME(vector<AD<AD<AD<double> > > > x){	\
  return NAME##checkpoint_namespace::eval(x);				\
}
