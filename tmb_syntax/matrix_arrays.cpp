// Shows use of vector, matrix and array operations.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data objects
  DATA_INTEGER(i);		// Scalar integer

  // Parameter objects
  PARAMETER(p)

  // Objects of double type so that they can be printed
  vector<double> v1(2); 
  v1 << 9,11;
  std::cout << "v1 = \n" 	<< v1 				<< "\n";
  matrix<double> m1(2,2); 
  m1 << 1,2,3,4;
  std::cout << "m1 = \n" 	<< m1 				<< "\n";
  array<double> a1(2,2); 
  a1 << 11,12,13,14;
  std::cout << "a1 = \n" 	<< a1 				<< "\n";

  // Obtaining dimensions of objects
  vector<int> d1 = a1.dim;		
  std::cout << "Dimension of array a1= \n" 	<< d1 				<< "\n";
  vector<int> d2(2);			
  d2(0) = m1.rows();
  d2(1) = m1.cols();
  std::cout << "Dimension of matrix m1= \n" 	<< d2 				<< "\n";
  int d3 = v1.size();			
  std::cout << "Dimension of vector v1= " 	<< d3 				<< "\n";

  // Matrix multiplication versus elementwise operations
  std::cout << "Element-wise product of arrays = \n"	<< a1*a1 				<< "\n";
  std::cout << "Element-wise division of arrays = \n" 	<< a1/a1 				<< "\n";
  std::cout << "Element-wise product of matrices = \n"	<< m1.array()*m1.array() 		<< "\n";
  std::cout << "matrix product = \n" 			<< m1*m1 				<< "\n";
  std::cout << "matrix-vector product = \n" 		<< m1*v1 				<< "\n";
  std::cout << "Element-wise product of vectors = \n" 	<< v1*v1 				<< "\n";
  std::cout << "Element-wise division of vectors = \n" 	<< v1/v1 				<< "\n";
  std::cout << "Inner product of vectors = " 		<< (v1*v1).sum()			<< "\n";

  // Indexing objects
  std::cout << "Element (1,1) of matrix m1= " 		<< m1(1,1) 				<< "\n";
  std::cout << "2nd row of matrix m1= \n" 		<< m1.row(1) 				<< "\n";
  std::cout << "2nd col of matrix m1= \n" 		<< m1.col(1) 				<< "\n";
  std::cout << "Element (1,1) of array a1= " 		<< a1(1,1) 				<< "\n";
  std::cout << "2nd row of array a1= \n"		<< a1.transpose().col(1)		<< "\n";
  std::cout << "2nd col of array a1= \n" 		<< a1.col(1) 				<< "\n";
  std::cout << "2nd element of vector v1= "		<< v1(1) 				<< "\n";

  // Subsetting matrices and vectors
  std::cout << "First element of v1= "			<< v1.head(1) 				<< "\n";
  std::cout << "Last element of v1= "			<< v1.tail(1) 				<< "\n";
  std::cout << "Block of m1 consisting of m1(1,1)= "	<< m1.block(0,0,1,1)			<< "\n";
  
  

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
