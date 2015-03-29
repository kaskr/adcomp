/* Simple configuration - does not require maintainance on the R side. 
   Configuration variables can either be set from R with config() or
   from the user template after the include statement.
*/

static struct config_struct{
  /* Configuration variables.
     - Default values _must_be_ specified with SET(var,value) 
     - Can be either bool or integer.
  */
  struct {
    bool all;        /* Print tracing information */
    bool parallel;   /* Trace info from parallel for loops */
    bool optimize;
    bool atomic;     /* Trace construction of atomic functions */
  } trace;
  struct {
    bool instantly;  /* Always optimize just after tape creation */
    bool parallel;   /* Allow optimize in parallel (memory consuming) */
  } optimize;
  struct {
    bool parallel;   /* Parallel tape creation */
  } tape;
  struct {
    bool all;        /* Print debug info */
    bool getListElement;
  } debug;

  int cmd;
  SEXP envir;
  void set(const char* name, bool &var, bool default_value){
    // cmd=0: set defaults in this struct.
    // cmd=1: copy from this struct to R.
    // cmd=2: copy from R to this struct.
    if(cmd==0)var=default_value;
    if(cmd==1)defineVar(install(name),asSEXP(var),envir);
    if(cmd==2)var=INTEGER(findVar(install(name),envir))[0];
  }
#define SET(name,value)set(#name,name,value);
  void set(){
    SET(trace.all,true);
    SET(trace.parallel,true);
    SET(trace.optimize,true);
    SET(trace.atomic,true);
    SET(debug.all,false);
    SET(debug.getListElement,false);
    SET(optimize.instantly,true);
    SET(optimize.parallel,false);
    SET(tape.parallel,true);
  }
#undef SET
  config_struct(){
    cmd=0;
    set();
  };
} config;

extern "C"
{
  SEXP TMBconfig(SEXP envir, SEXP cmd){
    config.cmd=INTEGER(cmd)[0];
    config.envir=envir;
    config.set();
    return R_NilValue;
  }
}
