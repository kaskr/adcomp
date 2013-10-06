/* To be removed */
#define RCPPAD_DEBUG 0
#define RCPPAD_PRINT(x)std::cout << #x << ": " << x << "\n"; std::cout.flush();

/* Turn on debug for Eigen ? */
#ifdef RCPPAD_SAFEBOUNDS
#undef NDEBUG
#undef eigen_assert
void eigen_Rprintf(const char* x);
#define eigen_assert(x) if (!(x)) { eigen_Rprintf("RcppAD has received an error from Eigen. "); \
                                  eigen_Rprintf("The following condition was not met:\n");          \
                                  eigen_Rprintf(#x);                                                \
                                  eigen_Rprintf("\nPlease check your matrix-vector bounds etc., "); \
                                  eigen_Rprintf("or run your program through a debugger.\n");       \
				  abort();}
#else
#undef NDEBUG
#define NDEBUG 1
#endif
#include <Eigen/Dense>
#include <Eigen/Sparse>
/* R must come after Eigen because conflict with length macro on mac os x */
#include <R.h>
#include <Rdefines.h>
void eigen_Rprintf(const char* x){Rprintf(x);}
/* Always turn off debug for cppad */
#undef NDEBUG
#define NDEBUG 1
#include "cppad/cppad.hpp"
#undef NDEBUG
#include "my/my.cpp"
using my::matrix;
using my::vector;
using CppAD::AD;
using CppAD::ADFun;
#include "convert.hpp" // asSEXP, asMatrix, asVector
#include "config.hpp"
#include "dnorm.hpp"   // harmless
#include "lgamma.hpp"  // harmless
#include "Vectorize.hpp"
#include "start_parallel.hpp"
#include "convenience.hpp"

/* Memory manager:
   Count the number of external pointers alive.
   When total number is zero it is safe to dyn.unload
   the library.
*/
static struct memory_manager_struct{
  int counter;
  void RegisterCFinalizer(){
    counter++;
  }
  void CallCFinalizer(){
    counter--;
  }
  memory_manager_struct(){
    counter=0;
  }
} memory_manager;
extern "C"{
  SEXP get_number_of_external_pointers_alive(){
    SEXP ans;
    PROTECT(ans = NEW_INTEGER(1));
    INTEGER(ans)[0]=memory_manager.counter;
    UNPROTECT(1);
    return ans;
  }
}

#ifdef _OPENMP
bool _openmp=true;
#else
bool _openmp=false;
#endif

/* Utility: Compile time test for Type=double */
template<class Type>
struct isDouble{
  enum{value=false};
};
template<>
struct isDouble<double>{
  enum{value=true};
};

/* Macros to obtain data and parameters from R */
#define PARAMETER_MATRIX(name) matrix<Type> name(asMatrix<Type>(	\
        getListElement(objective_function::parameters,#name)));		\
        objective_function::fill(name,#name);
#define PARAMETER_VECTOR(name) vector<Type> name(objective_function::fillShape(asVector<Type>(objective_function::getShape(#name)),#name));
#define PARAMETER(name) Type name(objective_function::fillShape(asVector<Type>(objective_function::getShape(#name)),#name)[0]);
#define DATA_VECTOR(name) vector<Type> name(asVector<Type>(	\
        getListElement(objective_function::data,#name)));
#define DATA_MATRIX(name) matrix<Type> name(asMatrix<Type>(	\
        getListElement(objective_function::data,#name)));
#define DATA_SCALAR(name) Type name(asVector<Type>(		\
        getListElement(objective_function::data,#name))[0]);
#define DATA_INTEGER(name) int name(CppAD::Integer(asVector<Type>(	\
	getListElement(objective_function::data,#name))[0]));
#define DATA_FACTOR(name) vector<int> name(asVector<int>(	\
        getListElement(objective_function::data,#name)));
#define NLEVELS(name) LENGTH(getAttrib(getListElement(this->data,#name),install("levels")))
#define DATA_SPARSE_MATRIX(name) Eigen::SparseMatrix<Type> name(my::asSparseMatrix<Type>( \
        getListElement(objective_function::data,#name)));
// REPORT() construct new SEXP so never report in parallel!
#define REPORT(name) if(isDouble<Type>::value && this->current_parallel_region<0){          \
                        defineVar(install(#name),asSEXP(name),objective_function::report);}
#define ADREPORT(name) objective_function::reportvector.push(name,#name);
#define PARALLEL_REGION if(this->parallel_region())
#define DATA_ARRAY(name) my::array<Type> name(my::asArray<Type>(	\
        getListElement(objective_function::data,#name)));	
#define PARAMETER_ARRAY(name) my::array<Type> name(objective_function::fillShape(my::asArray<Type>(objective_function::getShape(#name)),#name));


// kasper: Not sure used anywhere
/* Get the hessian sparsity pattern of ADFun object pointer */
template<class Type>
matrix<int> HessianSparsityPattern(ADFun<Type> *pf){
  int n=pf->Domain();
  vector<bool> Px(n * n);
  for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < n; j++)
	Px[ i * n + j ] = false;
      Px[ i * n + i ] = true;
    }
  pf->ForSparseJac(n, Px);
  vector<bool> Py(1); Py[0]=true;
  return asMatrix(vector<int>(pf->RevSparseHes(n,Py)),n,n);
}


/* get the list element named str, or return NULL */ 
SEXP getListElement(SEXP list, const char *str) 
{
  if(config.debug.getListElement)std::cout << "getListElement: " << str << " ";
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol); 
  int i; 
  for (i = 0; i < length(list); i++) 
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) 
      {
	elmt = VECTOR_ELT(list, i); 
	break; 
      }
  if(config.debug.getListElement)std::cout << "Length: " << LENGTH(elmt) << " ";
  if(config.debug.getListElement)std::cout << "\n";
  return elmt; 
}


/* Do nothing if we are trying to tape non AD-types */
void Independent(vector<double> x){}

/* Used by ADREPORT */
template <class Type>
struct report_stack{
  vector<const char*> names;
  vector<int> namelength;
  vector<Type> result;
  void clear(){
    names.resize(0);
    namelength.resize(0);
    result.resize(0);
  }
  /* Make space for n new items of given name */
  void increase(int n, const char* name){
    names.conservativeResize(names.size()+1);
    names[names.size()-1]=name;
    namelength.conservativeResize(namelength.size()+1);
    namelength[namelength.size()-1]=n;
    result.conservativeResize(result.size()+n);
  }
  // push scalar
  void push(Type x, const char* name){
    increase(1,name);
    result[result.size()-1]=x;
  }
  // push vector, matrix or array
  template<class VectorType>
  void push(VectorType x, const char* name){
    int n=x.size();
    int oldsize=result.size();
    increase(n,name);
    for(int i=0;i<n;i++)result[oldsize+i]=x[i];
  }
  // Cast to vector
  operator vector<Type>(){
    return result;
  }
  /* Get names (with replicates) to R */
  SEXP reportnames()
  {
    int n=result.size();
    SEXP nam;
    PROTECT(nam=allocVector(STRSXP,n));
    int k=0;
    for(int i=0;i<names.size();i++){
      for(int j=0;j<namelength[i];j++){
	SET_STRING_ELT(nam,k,mkChar(names[i]));
	k++;
      }
    }
    UNPROTECT(1);
    return nam;
  }
};

template <class Type>
class objective_function
{
// private:
public:
  SEXP data;
  SEXP parameters;
  SEXP report;
  
  int index;
  vector<Type> theta;
  vector<const char*> thetanames;
  report_stack<Type> reportvector; //Used by "ADREPORT"
  bool reversefill; // used to find the parameter order in user template (not anymore - use pushParname instead)
  vector<const char*> parnames; // One name for each PARAMETER_ in user template
  void pushParname(const char* x){
    parnames.conservativeResize(parnames.size()+1);
    parnames[parnames.size()-1]=x;
  }

  /* ================== For parallel Hessian computation
     Need three different parallel evaluation modes:
     (1) *Parallel mode* where a parallel region is evaluated iff 
         current_parallel_region == selected_parallel_region
     (2) *Serial mode* where all parallel region tests are evaluated 
         to TRUE so that "PARALLEL_REGION" tests are effectively removed.
	 A negative value of "current_parallel_region" or "selected_parallel_region" 
	 is used to select this mode (the default).
     (3) *Count region mode* where statements inside "PARALLEL_REGION{...}"
         are *ignored* and "current_parallel_region" is increased by one each
	 time a parallel region is visited.
     NOTE: The macro "PARALLEL_REGION" is supposed to be defined as
           #define PARALLEL_REGION if(this->parallel_region())
	   where the function "parallel_region" does the book keeping.
   */
  bool parallel_ignore_statements;
  int current_parallel_region;       /* Identifier of a code-fragment of user template */
  int selected_parallel_region;      /* Consider _this_ code-fragment */
  int max_parallel_regions;          /* Max number of parallel region identifiers,
				        e.g. max_parallel_regions=omp_get_max_threads(); 
				        probably best in most cases. */
  bool parallel_region(){            /* Is this the selected parallel region ? */
    bool ans;
    if(current_parallel_region<0 || selected_parallel_region<0)return true; /* Serial mode */
    ans = (selected_parallel_region==current_parallel_region) && (!parallel_ignore_statements);
    current_parallel_region++;
    if(max_parallel_regions>0)current_parallel_region=current_parallel_region % max_parallel_regions;
    return ans;
  }
  int count_parallel_regions(){
    current_parallel_region=0;       /* reset counter */
    selected_parallel_region=0;
    parallel_ignore_statements=true; /* Do not evaluate stuff inside PARALLEL_REGION{...} */
    this->operator()();              /* Run through users code */
    if(max_parallel_regions>0)return max_parallel_regions;
    else
    return current_parallel_region;
  }
  void set_parallel_region(int i){   /* Select parallel region (from within openmp loop) */
    current_parallel_region=0;
    selected_parallel_region=i;
    parallel_ignore_statements=false;
  }


  /* data_ and parameters_ are R-lists containing R-vectors or R-matrices.
     report_ is an R-environment.  
     The elements of the vector "unlist(parameters_)" are filled into "theta"
     which contains the default parameter-values. This happens during the
     *construction* of the objective_function object.
     The user defined template "objective_function::operator()" is called
     from "MakeADFunObject" which tapes the operations and creates the final
     ADFun-object.
  */
  objective_function(SEXP data_, SEXP parameters_, SEXP report_)
  {
    report=report_;
    data=data_;
    parameters=parameters_;
    /* Fill theta with the default parameters. 
       Pass R-matrices column major. */
    theta.resize(nparms(parameters_));
    index=0;
    int counter=0;
    SEXP obj=parameters_;
    for(int i=0;i<length(obj);i++){
      for(int j=0;j<length(VECTOR_ELT(obj,i));j++)
	{
	  theta[counter++]=Type(REAL(VECTOR_ELT(obj,i))[j]);
	}
    }
    thetanames.resize(theta.size());
    current_parallel_region=-1;
    selected_parallel_region=-1;
    max_parallel_regions=-1;
    reversefill=false;
  }
  SEXP defaultpar()
  {
    int n=theta.size();
    SEXP res;
    SEXP nam;
    PROTECT(res=allocVector(REALSXP,n));
    PROTECT(nam=allocVector(STRSXP,n));
    for(int i=0;i<n;i++){
      //REAL(res)[i]=CppAD::Value(theta[i]);
      REAL(res)[i]=value(theta[i]);
      SET_STRING_ELT(nam,i,mkChar(thetanames[i]));
    }
    setAttrib(res,R_NamesSymbol,nam);
    UNPROTECT(2);
    return res;
  }
  SEXP parNames()
  {
    int n=parnames.size();
    SEXP nam;
    PROTECT(nam=allocVector(STRSXP,n));
    for(int i=0;i<n;i++){
      SET_STRING_ELT(nam,i,mkChar(parnames[i]));
    }
    UNPROTECT(1);
    return nam;
  }
  
  /* FIXME: "Value" should be "var2par" I guess 
     kasper: Why not use asDouble defined previously? */
  double value(double x){return x;}
  double value(AD<double> x){return CppAD::Value(x);}
  double value(AD<AD<double> > x){return CppAD::Value(CppAD::Value(x));}
  
  int nparms(SEXP obj)
  {
    int count=0;
    for(int i=0;i<length(obj);i++){
      if(!isReal(VECTOR_ELT(obj,i)))error("PARAMETER COMPONENT NOT A VECTOR!");
      count+=length(VECTOR_ELT(obj,i));
    }
    return count;
  }

  /* The "fill functions" are all used to populate parameter vectors,
     arrays, matrices etc with the values of the parameter vector theta. */
  void fill(vector<Type> &x, const char *nam)
  {
    pushParname(nam);
    for(int i=0;i<x.size();i++){
      thetanames[index]=nam;
      if(reversefill)theta[index++]=x[i];else x[i]=theta[index++];
    }
  }
  void fill(matrix<Type> &x, const char *nam)
  {
    pushParname(nam);
    for(int j=0;j<x.cols();j++){
      for(int i=0;i<x.rows();i++){
	thetanames[index]=nam;
	if(reversefill)theta[index++]=x(i,j);else x(i,j)=theta[index++];
      }
    }
  }
  template<class ArrayType>
  void fill(ArrayType &x, const char *nam)
  {
    pushParname(nam);
    for(int i=0;i<x.size();i++){
	thetanames[index]=nam;
	if(reversefill)theta[index++]=x[i];else x[i]=theta[index++];
    }
  }

  /* Experiment: new map feature - currently arrays only */
  template<class ArrayType>
  void fillmap(ArrayType &x, const char *nam)
  {
    pushParname(nam);
    SEXP elm=getListElement(parameters,nam);
    int* map=INTEGER(getAttrib(elm,install("map")));
    int  nlevels=INTEGER(getAttrib(elm,install("nlevels")))[0];
    for(int i=0;i<x.size();i++){
      if(map[i]>=0){
	thetanames[index+map[i]]=nam;
	if(reversefill)theta[index+map[i]]=x[i];else x[i]=theta[index+map[i]];
      }
    }
    index+=nlevels;
  }
  // Auto detect whether we are in "map-mode"
  SEXP getShape(const char *nam){
    SEXP elm=getListElement(parameters,nam);
    SEXP shape=getAttrib(elm,install("shape"));
    if(shape==R_NilValue)return elm;
    else return shape;
  }
  template<class ArrayType>
  //ArrayType fillShape(ArrayType &x, const char *nam){
  ArrayType fillShape(ArrayType x, const char *nam){
    SEXP elm=getListElement(parameters,nam);
    SEXP shape=getAttrib(elm,install("shape"));
    if(shape==R_NilValue)fill(x,nam);
    else fillmap(x,nam);
    return x;
  }

  void fill(Type &x, char const *nam)
  {
    pushParname(nam);
    thetanames[index]=nam;
    if(reversefill)theta[index++]=x;else x=theta[index++];
  }
   
   Type operator() ();
};

/* Experiment: Help manage the parallel accumulation.
   Usage:
   parallel_accumulator<Type> res(this);
   For safety we avoid adding an operator=() !
 */
template<class Type>
struct parallel_accumulator{
  Type result;
  objective_function<Type>* obj;
  parallel_accumulator(objective_function<Type>* obj_){
    result=Type(0);
    obj=obj_;
#ifdef _OPENMP
#include <omp.h>
    obj->max_parallel_regions=omp_get_max_threads();
#endif
  }
  inline void operator+=(Type x){
    if(obj->parallel_region())result+=x;
  }
  inline void operator-=(Type x){
    if(obj->parallel_region())result-=x;
  }
  operator Type(){
    return result;
  }
};

/* Template to evaluate an ADFun object from R:
   Template argument can be
   - ADFun or
   - parallelADFun
*/
template<class ADFunType>
SEXP EvalADFunObjectTemplate(SEXP f, SEXP theta, SEXP control)
{
  if(!isNewList(control))error("'control' must be a list");
  ADFunType* pf;
  pf=(ADFunType*)R_ExternalPtrAddr(f);
  PROTECT(theta=AS_NUMERIC(theta));
  int n=pf->Domain();
  int m=pf->Range();
  if(LENGTH(theta)!=n)error("Wrong parameter length.");
  //R-index -> C-index
  int rangecomponent=INTEGER(getListElement(control,"rangecomponent"))[0]-1;
  if(!((0<=rangecomponent)&(rangecomponent<=m-1)))
    error("Wrong range component.");
  int order = INTEGER(getListElement(control,"order"))[0];
  if((order!=0) & (order!=1) & (order!=2) & (order!=3))
    error("order can be 0, 1, 2 or 3");
  int sparsitypattern=INTEGER(getListElement(control,"sparsitypattern"))[0];
  SEXP hessiancols; // Hessian columns
  PROTECT(hessiancols=getListElement(control,"hessiancols"));
  int ncols=length(hessiancols);
  SEXP hessianrows; // Hessian rows
  PROTECT(hessianrows=getListElement(control,"hessianrows"));
  int nrows=length(hessianrows);
  if((nrows>0)&(nrows!=ncols))error("hessianrows and hessianrows must have same length");
  vector<int> cols(ncols);    
  vector<int> cols0(ncols);
  vector<int> rows(nrows);
  if(ncols>0){
    for(int i=0;i<ncols;i++){
      cols[i]=INTEGER(hessiancols)[i]-1; //R-index -> C-index
      cols0[i]=0;
      if(nrows>0)rows[i]=INTEGER(hessianrows)[i]-1; //R-index -> C-index
    }
  }
  vector<double> x(n);
  for(int i=0;i<n;i++)x[i]=REAL(theta)[i];
  SEXP res=R_NilValue;
  SEXP rangeweight=getListElement(control,"rangeweight");
  if(rangeweight!=R_NilValue){
    if(LENGTH(rangeweight)!=m)error("rangeweight must have length equal to range dimension");
    pf->Forward(0,x);
    res=asSEXP(pf->Reverse(1,asVector<double>(rangeweight)));
    UNPROTECT(3);
    return res;
  }
  if(order==3){
    vector<double> w(1);
    w[0]=1;
    if((nrows!=1) | (ncols!=1))error("For 3rd order derivatives a single hessian coordinate must be specified.");
    pf->ForTwo(x,rows,cols); /* Compute forward directions */
    PROTECT(res=asSEXP(asMatrix(pf->Reverse(3,w),n,3)));
  }
  if(order==0){
    PROTECT(res=asSEXP(pf->Forward(rangecomponent,x)));
    SEXP rangenames=getAttrib(f,install("range.names"));
    if(LENGTH(res)==LENGTH(rangenames)){
      setAttrib(res,R_NamesSymbol,rangenames);
    }
  }
  if(order==1)PROTECT(res=asSEXP(asMatrix(pf->Jacobian(x),m,n)));
  //if(order==2)res=asSEXP(pf->Hessian(x,0),1);
  if(order==2){
    if(ncols==0){
      if(sparsitypattern){
	PROTECT(res=asSEXP(HessianSparsityPattern(pf)));  
      } else {
	PROTECT(res=asSEXP(asMatrix(pf->Hessian(x,rangecomponent),n,n)));
      }
    }
    else if (nrows==0){
      /* Fixme: the cols0 argument should be user changeable */
      PROTECT(res=asSEXP(asMatrix(pf->RevTwo(x,cols0,cols),n,ncols)));
    }
    else PROTECT(res=asSEXP(asMatrix(pf->ForTwo(x,rows,cols),m,ncols)));
  }
  UNPROTECT(4);
  return res;
}

/* How to garbage collect an ADFun or parallelADFun object pointer */
template <class ADFunType>
void finalize(SEXP x)
{
  ADFunType* ptr=(ADFunType*)R_ExternalPtrAddr(x);
  if(ptr!=NULL)delete ptr;
  memory_manager.CallCFinalizer();
}


ADFun<double>* MakeADFunObject(SEXP data, SEXP parameters,
			       SEXP report, SEXP control, int parallel_region=-1,
			       SEXP &info=R_NilValue)
{
  int returnReport = INTEGER(getListElement(control,"report"))[0];
  /* Create objective_function "dummy"-object */
  objective_function< AD<double> > F(data,parameters,report);
  F.set_parallel_region(parallel_region);
  /* Create ADFun pointer.
     We have the option to tape either the value returned by the
     objective_function template or the vector reported using the
     macro "ADREPORT" */
  Independent(F.theta);  // In both cases theta is the independent variable
  ADFun< double >* pf;
  if(!returnReport){ // Default case: no ad report - parallel run allowed
    vector< AD<double> > y(1);
    y[0]=F();
    pf = new ADFun< double >(F.theta,y);
  } else { // ad report case
    F(); // Run through user template (modifies reportvector)
    pf = new ADFun< double >(F.theta,F.reportvector.result);
    info=F.reportvector.reportnames(); // parallel run *not* allowed
  }
  return pf;
}


extern "C"
{

  /* How to garbage collect an ADFun object pointer */
  void finalizeADFun(SEXP x)
  {
    ADFun<double>* ptr=(ADFun<double>*)R_ExternalPtrAddr(x);
    if(ptr!=NULL)delete ptr;
    memory_manager.CallCFinalizer();
  }
  void finalizeparallelADFun(SEXP x)
  {
    parallelADFun<double>* ptr=(parallelADFun<double>*)R_ExternalPtrAddr(x);
    if(ptr!=NULL)delete ptr;
    memory_manager.CallCFinalizer();
  }

  /* Construct ADFun object */
  SEXP MakeADFunObject(SEXP data, SEXP parameters,
		       SEXP report, SEXP control)
  {
    /* Some type checking */
    if(!isNewList(data))error("'data' must be a list");
    if(!isNewList(parameters))error("'parameters' must be a list");
    if(!isEnvironment(report))error("'report' must be an environment");
    if(!isNewList(control))error("'control' must be a list");
    int returnReport = INTEGER(getListElement(control,"report"))[0];

    /* Get the default parameter vector (tiny overhead) */
    objective_function< AD<double> > F(data,parameters,report);
    F();
    SEXP par;
    PROTECT(par=F.defaultpar());
    SEXP res;
    SEXP info;
    PROTECT(info=R_NilValue); // Important

    if(_openmp && !returnReport){ // Parallel mode
#ifdef _OPENMP
      std::cout << "Count num parallel regions\n";
      objective_function< double > FF(data,parameters,report);
      int n=FF.count_parallel_regions(); // overhead: This could be taken from F instead
      std::cout << n << " regions found.\n";
      start_parallel(); /* Start threads */
      vector< ADFun<double>* > pfvec(n);
#pragma omp parallel for if (config.tape.parallel)
      for(int i=0;i<n;i++){
	pfvec[i]=MakeADFunObject(data, parameters, report, control, i, info);
	if(config.optimize.instantly)pfvec[i]->optimize();
      }
      parallelADFun<double>* ppf=new parallelADFun<double>(pfvec);
      /* Convert parallel ADFun pointer to R_ExternalPtr */
      PROTECT(res=R_MakeExternalPtr((void*) ppf,mkChar("parallelADFun"),R_NilValue));
      R_RegisterCFinalizer(res,finalizeparallelADFun);
      memory_manager.RegisterCFinalizer();
#endif
    } else { // Serial mode
      /* Actual work: tape creation */
      ADFun<double>* pf=MakeADFunObject(data, parameters, report, control, -1, info);
      /* Convert ADFun pointer to R_ExternalPtr */
      PROTECT(res=R_MakeExternalPtr((void*) pf,mkChar("ADFun"),R_NilValue));
      setAttrib(res,install("range.names"),info);
      R_RegisterCFinalizer(res,finalizeADFun);
      memory_manager.RegisterCFinalizer();
    }

    /* Return list of external pointer and default-parameter */
    SEXP ans;
    PROTECT(ans=allocVector(VECSXP,2));
    SET_VECTOR_ELT(ans,0,res);
    SET_VECTOR_ELT(ans,1,par);
    UNPROTECT(4);

    return ans;
  }
  
  SEXP InfoADFunObject(SEXP f)
  {
    ADFun<double>* pf;
    pf=(ADFun<double>*)R_ExternalPtrAddr(f);
    SEXP ans,names;
    PROTECT(ans=allocVector(VECSXP,4));
    PROTECT(names=allocVector(STRSXP,4));
    SET_VECTOR_ELT(ans,0,asSEXP(int(pf->Domain())));
    SET_STRING_ELT(names,0,mkChar("Domain"));
    SET_VECTOR_ELT(ans,1,asSEXP(int(pf->Range())));
    SET_STRING_ELT(names,1,mkChar("Range"));
    SET_VECTOR_ELT(ans,2,asSEXP(int(pf->use_VecAD())));
    SET_STRING_ELT(names,2,mkChar("use_VecAD"));
    SET_VECTOR_ELT(ans,3,asSEXP(int(pf->size_var())));
    SET_STRING_ELT(names,3,mkChar("size_var"));
    setAttrib(ans,R_NamesSymbol,names);
    UNPROTECT(2);
    return ans;
  }

  SEXP optimizeADFunObject(SEXP f)
  {
    SEXP tag=R_ExternalPtrTag(f);
    if(!strcmp(CHAR(tag), "ADFun")){
      ADFun<double>* pf;
      pf=(ADFun<double>*)R_ExternalPtrAddr(f);
      pf->optimize();
    }
    if(!strcmp(CHAR(tag), "parallelADFun")){
      parallelADFun<double>* pf;
      pf=(parallelADFun<double>*)R_ExternalPtrAddr(f);
      pf->optimize();      
    }
    return R_NilValue;
  }

  /* Get tag of external pointer */
  SEXP getTag(SEXP f){
    return R_ExternalPtrTag(f);
  }

  SEXP EvalADFunObject(SEXP f, SEXP theta, SEXP control)
  {
    SEXP tag=R_ExternalPtrTag(f);
    if(!strcmp(CHAR(tag), "ADFun"))
      return EvalADFunObjectTemplate<ADFun<double> >(f,theta,control);
    if(!strcmp(CHAR(tag), "parallelADFun"))
      return EvalADFunObjectTemplate<parallelADFun<double> >(f,theta,control);
    error("NOT A KNOWN FUNCTION POINTER");
  }
  
}

/* Double interface */
extern "C"
{

  /* How to garbage collect a DoubleFun object pointer */
  void finalizeDoubleFun(SEXP x)
  {
    objective_function<double>* ptr=(objective_function<double>*)R_ExternalPtrAddr(x);
    if(ptr!=NULL)delete ptr;
    memory_manager.CallCFinalizer();
  }
  
  SEXP MakeDoubleFunObject(SEXP data, SEXP parameters, SEXP report)
  {
    /* Some type checking */
    if(!isNewList(data))error("'data' must be a list");
    if(!isNewList(parameters))error("'parameters' must be a list");
    if(!isEnvironment(report))error("'report' must be an environment");
    
    /* Create DoubleFun pointer */
    objective_function<double>* pF = 
      new objective_function<double>(data,parameters,report);
    
    /* Convert DoubleFun pointer to R_ExternalPtr */
    SEXP res;
    PROTECT(res=R_MakeExternalPtr((void*) pF,mkChar("DoubleFun"),R_NilValue));
    R_RegisterCFinalizer(res,finalizeDoubleFun);
    memory_manager.RegisterCFinalizer();
    UNPROTECT(1);

    return res;
  }

  
  SEXP EvalDoubleFunObject(SEXP f, SEXP theta, SEXP control)
  {
    objective_function<double>* pf;
    pf=(objective_function<double>*)R_ExternalPtrAddr(f);
    PROTECT(theta=AS_NUMERIC(theta));
    int n=pf->theta.size();
    if(LENGTH(theta)!=n)error("Wrong parameter length.");

    vector<double> x(n);
    for(int i=0;i<n;i++)x[i]=REAL(theta)[i];
    pf->theta=x;
    
    /* Since we are actually evaluating objective_function::operator() (not
       an ADFun object) we should remember to initialize parameter-index. */
    pf->index=0;
    pf->parnames.resize(0); // To avoid mem leak.
    pf->reportvector.clear();
    SEXP res;
    res=asSEXP(pf->operator()());
    UNPROTECT(1);
    return res;
  }

  /* We spend a function evaluation on getting the parameter order (!) */
  SEXP getParameterOrder(SEXP data, SEXP parameters, SEXP report)
  {
    /* Some type checking */
    if(!isNewList(data))error("'data' must be a list");
    if(!isNewList(parameters))error("'parameters' must be a list");
    if(!isEnvironment(report))error("'report' must be an environment");
    objective_function<double> F(data,parameters,report);
    //F.reversefill=true;
    F(); // Run through user template
    //return(F.defaultpar());
    return F.parNames();
  }

}



extern "C"
{

  /* Tape the gradient using nested AD types */
  SEXP MakeADGradObject(SEXP data, SEXP parameters, SEXP report)
  {
    /* Some type checking */
    if(!isNewList(data))error("'data' must be a list");
    if(!isNewList(parameters))error("'parameters' must be a list");
    if(!isEnvironment(report))error("'report' must be an environment");

    /* Create ADFun pointer */
    objective_function< AD<AD<double> > > F(data,parameters,report);
    int n=F.theta.size();
    Independent(F.theta);
    vector< AD<AD<double> > > y(1);
    y[0]=F();
    ADFun<AD<double> > tmp(F.theta,y);

    vector<AD<double> > x(n);
    for(int i=0;i<n;i++)x[i]=CppAD::Value(F.theta[i]);
    vector<AD<double> > yy(n);
    Independent(x);
    yy=tmp.Jacobian(x);
    ADFun< double >* pf = new ADFun< double >(x,yy);
    
    /* Get the default parameter vector */
    SEXP par;
    PROTECT(par=F.defaultpar());

    /* Convert ADFun pointer to R_ExternalPtr */
    SEXP res;
    PROTECT(res=R_MakeExternalPtr((void*) pf,mkChar("ADFun"),R_NilValue));
    R_RegisterCFinalizer(res,finalizeADFun);
    memory_manager.RegisterCFinalizer();

    /* Return list */
    SEXP ans;
    PROTECT(ans=allocVector(VECSXP,2));
    SET_VECTOR_ELT(ans,0,res);
    SET_VECTOR_ELT(ans,1,par);
    UNPROTECT(3);

    return ans;
  }
}


/* ======================== EXPERIMENT: Tape sparse hessian */
/* Status: Very effecient to evaluate (after optimize) however tape creation takes too long */
extern "C"
{

  /* Tape the hessian[cbind(i,j)] using nested AD types */
  SEXP MakeADHessObject(SEXP data, SEXP parameters, SEXP report, SEXP hessianrows, SEXP hessiancols)
  {
    /* Some type checking */
    if(!isNewList(data))error("'data' must be a list");
    if(!isNewList(parameters))error("'parameters' must be a list");
    if(!isEnvironment(report))error("'report' must be an environment");
    
    /* TODO: Check i anf j*/
    int ncols=length(hessiancols);
    int nrows=length(hessianrows);
    if(nrows!=ncols)error("hessianrows and hessiancols must have same length");
    int m=ncols; /* Dimension of range vector */
    vector<int> cols(ncols);    
    vector<int> rows(nrows);
    for(int i=0;i<ncols;i++){
      cols[i]=INTEGER(hessiancols)[i]-1; //R-index -> C-index
      rows[i]=INTEGER(hessianrows)[i]-1; //R-index -> C-index
    }
    
    /* Create ADFun pointer */
    objective_function< AD<AD<double> > > F(data,parameters,report);
    int n=F.theta.size();
    Independent(F.theta);
    vector< AD<AD<double> > > y(1);
    y[0]=F();
    ADFun<AD<double> > tmp(F.theta,y);
    
    /* optimize it */
    tmp.optimize();

    vector<AD<double> > x(n);
    for(int i=0;i<n;i++)x[i]=CppAD::Value(F.theta[i]);
    vector<AD<double> > yy(m);
    Independent(x);
    //yy=tmp.Jacobian(x);
    yy=tmp.ForTwo(x,rows,cols);
    ADFun< double >* pf = new ADFun< double >(x,yy);
    
    /* Get the default parameter vector */
    SEXP par;
    PROTECT(par=F.defaultpar());

    /* Convert ADFun pointer to R_ExternalPtr */
    SEXP res;
    PROTECT(res=R_MakeExternalPtr((void*) pf,mkChar("ADFun"),R_NilValue));
    R_RegisterCFinalizer(res,finalizeADFun);
    memory_manager.RegisterCFinalizer();

    /* Return list */
    SEXP ans;
    PROTECT(ans=allocVector(VECSXP,2));
    SET_VECTOR_ELT(ans,0,res);
    SET_VECTOR_ELT(ans,1,par);
    UNPROTECT(3);

    return ans;
  }



}


/* ======================== EXPERIMENT: 
   Tape gradient on AD<AD<double>> and run optimize method. 
   Then tape the sparse hessian as the gradient of each component.
 */
  /* Tape the hessian[cbind(i,j)] using nested AD types 

     skip: integer vector of columns to skip from the hessian (will not change dimension
            - only treat h[:,skip] and h[skip,:] as zero). Negative subscripts are not allowed.
   */
sphess MakeADHessObject2(SEXP data, SEXP parameters, SEXP report, SEXP skip, int parallel_region=-1)
{
  /* Some type checking */
  if(!isNewList(data))error("'data' must be a list");
  if(!isNewList(parameters))error("'parameters' must be a list");
  if(!isEnvironment(report))error("'report' must be an environment");
  
  int m;  
  //int skip_=INTEGER(skip)[0];
  //#define KEEP(i)(i<n-skip_)
  //#define KEEP(i)(i>=skip_)
  /* Keep lower triangle of hessian - and only random effect part */
  //#define KEEP_COL(col)(col<n-skip_)
  //#define KEEP_ROW(row,col)(KEEP_COL(row)&(row>=col))
  /* Create ADFun pointer */
  objective_function< AD<AD<AD<double> > > > F(data,parameters,report);
  F.set_parallel_region(parallel_region);
  int n=F.theta.size();

  vector<bool> keepcol(n); // Scatter for fast lookup 
  for(int i=0;i<n;i++){keepcol[i]=true;}
  for(int i=0;i<LENGTH(skip);i++){
    keepcol[INTEGER(skip)[i]-1]=false; // skip is R-index !
  }
#define KEEP_COL(col)(keepcol[col])
#define KEEP_ROW(row,col)(KEEP_COL(row)&(row>=col))

  Independent(F.theta);
  vector< AD<AD<AD<double> > > > y(1);
  y[0]=F();
  ADFun<AD<AD<double> > > tmp(F.theta,y);

  /* Tape gradient R^n -> R^n */
  vector<AD<AD<double> > > xx(n);
  for(int i=0;i<n;i++)xx[i]=CppAD::Value(F.theta[i]);
  vector<AD<AD<double> > > yy(n);
  Independent(xx);
  yy=tmp.Jacobian(xx);
  ADFun<AD<double > > tmp2(xx,yy);
  
  /* Optimize tape */
  tmp2.optimize();  // ================== WARNING!!!
  
  /* ========================================================== NOT DONE YET */
  /* Tape hessian  */
  tmp2.my_init();
  //std::cout << tmp2.colpattern << "\n";
  m=0;
  //for(int i=0;i<tmp2.colpattern.size();i++)m+=tmp2.colpattern[i].size();
  int colisize;
  for(int i=0;i<int(tmp2.colpattern.size());i++){
    colisize=tmp2.colpattern[i].size();
    if(KEEP_COL(i)){
      for(int j=0;j<colisize;j++){
	m+=KEEP_ROW( tmp2.colpattern[i][j] , i);
      }
    }
  }
    
  // argument and result for reverse mode calculations
  vector<AD<double> > u(n);
  vector<AD<double> > v(n);
    
    
  // initialize all the components
  for(int i = 0; i < n; i++)
    v[i] = 0.0;
  
  vector<AD<double> > xxx(n);
  for(int i=0;i<n;i++)xxx[i]=CppAD::Value(CppAD::Value(F.theta[i]));
  vector<AD<double> > yyy(m);
  CppAD::vector<int>* icol;
  int k=0;
  
  Independent(xxx);
  // Take from Jacobian.hpp ...
  // point at which we are evaluating the Jacobian
  tmp2.Forward(0, xxx);
  for(int i = 0; i < n; i++){
    // RCPPAD_PRINT(i);
    /* ========== ORIGINAL
       v[i] = 1.0;
       u = tmp2.Reverse(1, v);
       v[i] = 0.0;
       yyy[i]=u[i]; // <--- Test: return the diagonal entry. 
       ========== ORIGINAL */
    if(KEEP_COL(i)){
      tmp2.myReverse(1, v, i /*range comp*/, u /*domain*/);
      icol=&tmp2.colpattern[i];
      for(int j=0;j<int(icol->size());j++){
	if(KEEP_ROW( icol->operator[](j), i ))yyy[k++]=u[icol->operator[](j)];
      }
    }
    
  }
  
  //yyy=tmp2.Jacobian(xxx);

  //    ADFun< double >* pf = new ADFun< double >(xxx,yyy);
  ADFun< double >* pf = new ADFun< double >;
  pf->Dependent(xxx,yyy);

  if(config.optimize.instantly){
    if(config.trace.optimize)std::cout << "Optimizing tape... ";
    pf->optimize();
    if(config.trace.optimize)std::cout << "Done\n";
  }
  /* ========================================================= */    
  
  /* Calculate row and col index vectors.
     The variable "k" now holds the number of non-zeros.
  */

  // colproj used to remap (i,j) pairs. Not used anymore.
  vector<int> colproj(n); /* column projection index vector */
  int cumsum=0;
  for(int i = 0; i < n; i++){colproj[i]=cumsum;cumsum+=KEEP_COL(i);}
  
  vector<int> rowindex(k);
  vector<int> colindex(k);
  k=0;
  for(int i = 0; i < n; i++){
    if(KEEP_COL(i)){
      icol=&tmp2.colpattern[i];
      for(int j=0;j<int(icol->size());j++){
	if(KEEP_ROW( icol->operator[](j), i )){
	  // rowindex[k]=colproj[icol->operator[](j)];
	  // colindex[k]=colproj[i];
	  rowindex[k]=icol->operator[](j);
	  colindex[k]=i;
	  k++;
	}
      }
    }
  }
  sphess ans(pf,rowindex,colindex);
  return ans;
}

// kasper: Move to new file e.g. "convert.hpp"
template <class ADFunType>
SEXP asSEXP(const sphess_t<ADFunType> &H, const char* tag)
{
    SEXP par;
    par=R_NilValue;
    /* Convert ADFun pointer to R_ExternalPtr */
    SEXP res;
    //PROTECT(res=R_MakeExternalPtr((void*) H.pf,mkChar("ADFun"),R_NilValue));
    PROTECT(res=R_MakeExternalPtr((void*) H.pf,mkChar(tag),R_NilValue));


    //R_RegisterCFinalizer(res,finalizeADFun);
    R_RegisterCFinalizer(res,finalize<ADFunType>);
    memory_manager.RegisterCFinalizer();

    /* Return list */
    SEXP ans;
    PROTECT(ans=allocVector(VECSXP,4));
    SET_VECTOR_ELT(ans,0,res);
    SET_VECTOR_ELT(ans,1,par);
    SET_VECTOR_ELT(ans,2,asSEXP(H.i));
    SET_VECTOR_ELT(ans,3,asSEXP(H.j));
    UNPROTECT(2);
    return ans;
}




extern "C"
{
#ifdef _OPENMP
  SEXP MakeADHessObject2(SEXP data, SEXP parameters, SEXP report, SEXP skip){

    std::cout << "Count num parallel regions\n";
    objective_function< double > F(data,parameters,report);
    int n=F.count_parallel_regions();
    std::cout << n << " regions found.\n";

    start_parallel(); /* Start threads */

    /* parallel test */
    vector<sphess*> Hvec(n);
#pragma omp parallel for if (config.tape.parallel)
    //for(int i=0;i<n;i++)Hvec[i]=&MakeADHessObject2(data, parameters, report, skip, i);
    for(int i=0;i<n;i++)Hvec[i]=new sphess(MakeADHessObject2(data, parameters, report, skip, i));

    //parallelADFun<double> tmp(Hvec);
    parallelADFun<double>* tmp=new parallelADFun<double>(Hvec);

    
    return asSEXP(tmp->convert(),"parallelADFun");

  }
#else
  SEXP MakeADHessObject2(SEXP data, SEXP parameters, SEXP report, SEXP skip){
    //MakeADHessObject2_parallel(data, parameters, report, skip);
    //sphess H=MakeADHessObject2(data, parameters, report, skip, 10000/* -1*/);
    sphess H=MakeADHessObject2(data, parameters, report, skip, -1);
    //RCPPAD_PRINT(H.pf->Range());
    /* Get the default parameter vector */
    return asSEXP(H,"ADFun");
  }
#endif
  
}
