// Shows use of vector, matrix and array operations.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data objects
  DATA_VECTOR(v1);		
  DATA_MATRIX(m1);		
  DATA_ARRAY(a1);		

  // Parameter objects
  PARAMETER(p)		// Not used in this example

  // Obtaining dimensions of objects
  REPORT(a1.dim);
  vector<int> d2(2);			
  d2(0) = m1.rows();
  d2(1) = m1.cols();
  REPORT(d2);
  int d3 = v1.size();			
  REPORT(d3);
  

  // Matrix multiplication versus elementwise multiplication (similar for addition, subtraction,...)
  matrix<Type> m1m1 = m1*m1;	// Matrix multiplication of matrices
  REPORT(m1m1);  
  array<Type> a1a1 = a1*a1;	    // Element-wise multiplication of arrays 
  REPORT(a1a1);  
  matrix<Type> m1m1_by_element = (m1.array()*m1.array()).matrix();	// Elementwise multiplication of matrices
  REPORT(m1m1_by_element);  
  array<Type> a1a1_matrix_mult(2,2);
  a1a1_matrix_mult = (a1.matrix()*a1.matrix()).array();	// Matrix multiplication of arrays
  REPORT(a1a1_matrix_mult);  

  // Matrix-vector multiplication
  REPORT(m1*v1);	// matrix-vector product (linear algebra style)
  Type v1_norm2 = (v1*v1).sum();	// Inner product v1*v1
  REPORT(v1_norm2);
  
  // Indexing objects
  m1(1,1); 		// Element (1,1) of matrix m1
  m1.row(1);    // 2nd row of matrix m1
  m1.col(1); 	// 2nd col of matrix m1
  a1(1,1);		// Element (1,1) of array a1
  a1.transpose().col(1); // 2nd row of array a1
  a1.col(1); 	//2nd col of array a1
  v1(1);		// 2nd element of vector v1

  // Subsetting matrices and vectors
  v1.head(1);		// First element of v1
  v1.tail(1); 		// Last element of v1
  m1.block(0,0,1,1);	// Block of m1 consisting of m1(1,1)

  // Subsetting arrays
  // See R help pages for "template", i.e. "help(template)" in R  
  

  // Generic matrix operations that we must ensure compiles
  m1.transpose();
  m1.diagonal(); 
  m1.asDiagonal();
  

  Type ans;
  return ans;

}

/** \file
\ingroup matrix_arrays Examples
\brief Shows the use of vector, matrix and array operations.

*/
