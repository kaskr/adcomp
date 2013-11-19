// \file matrix_arrays_file
// \brief Shows examples of vector, matrix and array operations
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data objects
  DATA_INTEGER(i);		// Scalar integer

  // Parameter objects
  PARAMETER(p)

  // Objects of double type so that they can be printed
  vector<double> v1(3); 
  v1.fill(2.0);
  matrix<double> m1(3,3); 
  m1.fill(3.0);
  array<double> a1(3,3); 
  a1.fill(4.0);

  //Indexing
  vector<int> d1 = a1.dim;		// d1 will have length 1,2,3... depending of number of dimensions of a1
  vector<int> d2(2);			// Will hold dimension of matrix m1
  d2(0) = m1.rows();
  d2(1) = m1.cols();
  int d3 = v1.size();			// Length of vector

  // Matrix multiplication versus elementwise operations
  std::cout << "Elementwise product of arrays = \n" 	<< a1*a1 				<< "\n";
  std::cout << "Elementwise product of matrices = \n"	<< m1.array()*m1.array() 		<< "\n";
  std::cout << "matrix product = \n" 			<< m1*m1 				<< "\n";
  std::cout << "matrix-vector product = \n" 		<< m1*v1 				<< "\n";
  std::cout << "Elementwise product of vectors = \n" 	<< v1*v1 				<< "\n";
  std::cout << "Elementwise product of vectors = \n" 	<< v1.matrix().transpose()*v1.matrix() 	<< "\n";


  Type ans;
  return ans;

}
