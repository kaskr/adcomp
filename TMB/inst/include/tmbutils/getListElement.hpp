/* Helpers, to check that data and parameters are of the right types.
   "RObjectTester" denotes the type of a pointer to a test function.
   Examples of test functions are "isMatrix", "Rf_isArray", "isNumeric",
   etc (see Rinternals.h).
*/
typedef Rboolean (*RObjectTester)(SEXP);
#ifdef WITH_LIBTMB
void RObjectTestExpectedType(SEXP x, RObjectTester expectedtype, const char *nam);
Rboolean isValidSparseMatrix(SEXP x);
Rboolean isNumericScalar(SEXP x);
#else
void RObjectTestExpectedType(SEXP x, RObjectTester expectedtype, const char *nam){
  if(expectedtype != NULL){
    if(!expectedtype(x)){
      if(Rf_isNull(x)){
	Rf_warning("Expected object. Got NULL.");
      }
      Rf_error("Error when reading the variable: '%s'. Please check data and parameters.",nam);
    }
  }
}
Rboolean isValidSparseMatrix(SEXP x){
  if(!Rf_inherits(x,"dgTMatrix"))Rf_warning("Expected sparse matrix of class 'dgTMatrix'.");
  return Rf_inherits(x,"dgTMatrix");
}
Rboolean isNumericScalar(SEXP x){
  if(LENGTH(x)!=1){
    Rf_warning("Expected scalar. Got length=%i",LENGTH(x));
    return FALSE;
  }
  return Rf_isNumeric(x);
}
#endif

/** \internal \brief Get list element named "str", or return NULL */
#ifdef WITH_LIBTMB
SEXP getListElement(SEXP list, const char *str, RObjectTester expectedtype=NULL);
int  getListInteger(SEXP list, const char *str, int default_value = 0);
#else
SEXP getListElement(SEXP list, const char *str, RObjectTester expectedtype=NULL)
{
  if(config.debug.getListElement)std::cout << "getListElement: " << str << " ";
  SEXP elmt = R_NilValue, names = Rf_getAttrib(list, R_NamesSymbol);
  int i; 
  for (i = 0; i < Rf_length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) 
      {
	elmt = VECTOR_ELT(list, i); 
	break; 
      }
  if(config.debug.getListElement)std::cout << "Length: " << LENGTH(elmt) << " ";
  if(config.debug.getListElement)std::cout << "\n";
  RObjectTestExpectedType(elmt, expectedtype, str);
  return elmt; 
}
int getListInteger(SEXP list, const char *str, int default_value = 0) {
  SEXP tmp = getListElement(list, str);
  if ( tmp == R_NilValue ) {
    Rf_warning("Missing integer variable '%s'. Using default: %d. (Perhaps you are using a model object created with an old TMB version?)", str, default_value);
    return default_value;
  }
  return INTEGER(tmp)[0];
}
#endif
