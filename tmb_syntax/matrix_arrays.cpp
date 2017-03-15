// Overview of vector, matrix and array operations in TMB.
// Use of the file: 
// - Browse through and find the construction you need.
// - I you want to inspect the result, use REPORT().
//   Example: if you uncomment the line "REPORT(v1_size);"
//            the variable v1_size will be returned and printed in R 

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(v1);
  DATA_VECTOR(v2);
  DATA_VECTOR(v3);
  DATA_MATRIX(m1);
  DATA_MATRIX(m2);
  DATA_MATRIX(m3);
  DATA_MATRIX(m4);
  DATA_ARRAY(a1);
  DATA_ARRAY(a2);

  PARAMETER(p);      
  
  // VECTOR Operations =====================================================================
  
  int v1_size = v1.size();              // Length of v1			
  //REPORT(v1_size);
  
  
  int v1_1 = 1;
  
  vector<Type> v2_head = v2.head(2);    // First n cells of v2 where n = 2
             
  vector<Type> v1_tail = v1.tail(v1_1); // Last n cells of v1 where n = v1_1

  vector<Type> v3_segment = v3.segment(2,3); // start at 3rd cell of v3 & take next 3 cells

  Type v1_sum = v1.sum();               // sum of all cells in v1
 
  Type sum_v1 = sum(v1);                // another way to sum all cells in v1 (often used to sum vector of likelihoods) 
 
  Type v1_prod = v1.prod();             // product of all cells in v1

  vector<Type> v1_plus_v1 = v1 + v1;    // cell by cell addition of two vectors (+, -, *, / available)
 
  vector<Type> v2_times_v2 = v2*v2;     // cell by cell multiplication of two vectors 

  Type v1_times_v1_sum = (v1*v1).sum(); // inner product of two vectors (sum 1xn * nX1 = 1X1)

  vector<Type> exp_v1 = exp(v1);        // exp of v1 (others available e.g. log)
 
  // Matrices ==============================================================================
  
  int m1_size = m1.size();              // number of cells in matrix (not dimension of matrix)
 
  vector<int> m1_dim(2);                // .dim does not work on matrices; note int vs Type
  m1_dim(0)            = m1.rows();
  m1_dim(1)            = m1.cols();
  
  Type m1_11 = m1(v1_1,0);              // cell (v1_1, 0) of m1
 
  int ii = CppAD::Integer(m1(0,0));     // int ii = m1(0,0); extraction of integers from matrices & arrays requires CppAD
  int jj = CppAD::Integer(m1(0,1));
  
  Type m2_iijj = m2(ii,jj);
  
  vector<Type> m1_row1 = m1.row(1);    // 2nd row of m1
  
  vector<Type> m1_col1 = m1.col(v1_1); // 2nd col of m1

  m1.col(1) = v1;                      // assign v1 to 2nd col of m1 (affects m1 based operations below!)
  
  Type m1_sum = m1.sum();              // sum of all cells in m1

  Type sum_m1 = sum(m1);               // another way to sum cells in m1         
 
  Type m1_prod = m1.prod();

  Type m1_mean = m1.mean();
 
  Type m1_mincoeff = m1.minCoeff();    // minimum value in m1

  Type m1_maxcoeff = m1.maxCoeff();    // maximum value in m1

  // Need to be careful using min & max operations on parameter dependent objects (see Wiki's Things-you-should-NOT-do-in-TMB)
  
  Type m2_trace = m2.trace();          // sum of m2 diagonal starting at (0,0)
  
  matrix<Type> m1_trans = m1.transpose(); // Transpose of m1

  vector<Type> m4_diag = m4.diagonal();   // Diagional of m1
  
  matrix<Type> m1_plus_m1_by_cell = m1 + m1;                             // cellwise + or - of matrices

  matrix<Type> m2_times_m2_by_cell = m2.array()*m2.array()*m2.array();	 // cellwise * or / of matrices                    

  // can chain commands eg matrix.colwise().sum().maxCoeff()
 
  vector<Type> m4_colwise_sum = m4.colwise().sum(); // m4 col sums

  vector<Type> m4_rowwise = m4.rowwise().sum();     // m4 row sums
 
  vector<Type> m4_colmean = m4.colwise().mean();    // m4 col mean
 
  // Math operators
  
  matrix<Type> m2_exp = exp(m2.array());  // exp(m1) does not work; must use array designation

  matrix<Type> m2_log = log(m2.array());  // log(m1) does not work; must use array designation

  matrix<Type> exp_m2_plus_m2_by_cell = exp(m2.array()) + m2.array();	  // when using exp, exp(m1.array()) + m2 doesn't work
        
  matrix<Type> exp_m2_times_m2_by_cell = exp(m2.array()) * m2.array();	// when using exp, exp(m1.array()) * m2 doesn't work
     
  // when chaining commands, need to identify with parentheses matrix being summarized
  
  matrix<Type> exp_m2_times_m2_sumed_by_row = (exp(m2.array())*m2.array()).matrix().rowwise().sum();	

  // combining vectors, matrices & arrays 
  
  matrix<Type> m1_times_v11 = m1 * v1(1);            // matrix * 2nd cell of v1

  matrix<Type> a1_plus_m1 = a1.matrix() + m1;        // a1 + m1.array() does not work
 
  // Matrix Algebra

  vector<Type> m1_times_v1 = m1*v1;                  // matrix-vector product 2x2 * 2x1 = 2x1 matrix
     
  matrix<Type> m1_times_m1 = m1*m1;	                 // matrix-matrix product 2X2 * 2x2 = 2x2 matrix
  
  matrix<Type> v2_matrix = v2.matrix();              // creates nX1 matrix from vector

  matrix <Type> v3_trans = v3.matrix().transpose();  // creates 1Xn matrix from vector

  Type exp_m2_1_sumed = exp(m2.row(1).array()).matrix() * v3.matrix();	
  
  matrix<Type> v1_v2_outer_prod = v1.matrix()*v2.matrix().transpose();   // Outer product (Eigen does not have special function)
  //REPORT(v1_v2_outer_prod);                                              // 2x1 * 1x4 = 2x4 matrix
  
  matrix<Type> exp_v1_v2_outer_prod = exp(v1).matrix()*exp(v2).matrix().transpose();   
  
 // extracting parts of matrices
 
  matrix<Type> m1_block = m1.block(0,0,1,1);                       // m1 block starting at m1(0,0) and taking 1 row & 1 col
       
  m2.block(1,1,1,2) = m1.block(0,0,1,2);                           // m1 block starting at m1(0,0) & taking 1 row & 2 cols  & assigning to m2 starting at m2(1,1) 
 
  m2.block(0,0,1,5) = v3_trans;                                    // assign v3 to first row of m2

  m2.block(1,0,1,3) = v3_trans.block(0,0,1,3) * Type(10.0);        // assign 10xv3 to first row of m2

  m2.block(1,1 ,1,m1_dim(1)) = m1.block(0,0,1,m1_dim(1));           // m1 block starting at m1(0,0) & taking 1 row & m1_dim(1) cols & assigning to m2 starting at m2(1,1) 
  
  vector<Type> m4_block_mean = m4.block(0,0,3,4).colwise().mean();  // col mean of block
 
// ARRAYS =====================================================================================
  
  int a1_size = a1.size();                                     // number of cells in array (not dimension)
   
  vector<int> a1_dim = a1.dim;                                 // dimension of a1
  
  int a1_rows = a1.rows();
  int a1_cols = a1.cols();
  
  Type a1_11 = a1(0,0);

  a1.col(1) = v1;                                               // assign v1 to 2nd col of a1
 
  Type a1_sum = a1.sum();                                       // sum of all cells in a1

  Type sum_a1 = sum(a1);                                        // another way to sum all cells in a1

  Type a1_mean = a1.mean();                                     // mean of all cells in a1
 
  array<Type> a2_transpose = a2.transpose();                    // transpose array a2
 
  array<Type> a2_rotate = a2.rotate(1);                         // rotate a2 n = one dimension (i.e. [7,5] to [5,7])

  vector<Type> a2_colwise_sum = a2.matrix().colwise().sum();    // a2.colwise().sum() adds all cells; need to convert to matrix & then do colwise

  vector<Type> a1_transpose_col1 = a1.transpose().col(1);       // 2nd row of array a1
      
  // cellwise * & / of arrays
  
  array<Type> a2_times_a2_times_a2_by_cell = a2*a2*a2;       

  // TMB loses dim attribute of array when using exp or log; need to explicitly define dimension attribute of array
  
  array<Type> a2_exp(a2.dim);                
  a2_exp = exp(a2);
  
  // array<Type> v1xv2_arr_outer_prod = v1.array()*v2.array().transpose();   doesn't work

  array<Type> a1a1_matrix_mult(2,2);
  a1a1_matrix_mult = (a1.matrix()*a1.matrix()).array();  // Matrix multiplication of arrays

  Type ans;
  return ans;
}
