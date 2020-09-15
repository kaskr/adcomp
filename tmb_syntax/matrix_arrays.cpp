// Overview of vector, matrix and array operations in TMB.

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

  // Vectors ========================================================

  int v1_size = v1.size();              // Length of v1
  REPORT(v1_size);

  vector<Type> v2_head = v2.head(2);    // First 2 elements of v2
  REPORT(v2_head);                      //  v2.tail(2): last 2 elem.

  vector<Type> v3_segment = v3.segment(2,3); //segment of 3 elements,
  REPORT(v3_segment);                     // starting at 3rd element

  vector<Type> v_concat_v1_v2(v1.size()+v2.size()); // Result vector must be allocated in advance!
  v_concat_v1_v2 << v1,v2;              // Syntax for concatination of v1 and v2
  REPORT(v_concat_v1_v2);

  Type v1_sum = v1.sum();               // sum of all cells in v1
  REPORT(v1_sum);

  Type sum_v1 = sum(v1);                // alternative summation
  REPORT(sum_v1);

  Type v1_prod = v1.prod();             // product of all cells in v1
  REPORT(v1_prod);

  vector<Type> v1_plus_v1 = v1 + v1;    // elementwise addition
  REPORT(v1_plus_v1);                   //   (-, *, / similar)

  Type v1_times_v1_sum = (v1*v1).sum(); // inner product
  REPORT(v1_times_v1_sum);

  vector<Type> exp_v1 = exp(v1);        // exp of v1 (also log)
  REPORT(exp_v1);

  // min() and max(). NB! Non-differentiable functions
  Type v1_mincoeff = min(v1);    // minimum value in v1
  REPORT(v1_mincoeff);

  matrix<Type> v2_matrix = v2.matrix(); //vector-to-matrix conversion
  REPORT(v2_matrix);

  // Matrices =======================================================

  int m1_size = m1.size();    // number of elements (not dimension)
  REPORT(m1_size);

  Type m1_11 = m1(1,0);       // Element (1, 0) of m1
  REPORT(m1_11);

  vector<Type> m1_row1 = m1.row(0);    // 1st row of m1
  REPORT(m1_row1);

  vector<Type> m1_col1 = m1.col(1);    // 2dn col of m1
  REPORT(m1_col1);

  m1.col(1) = v1;                      // assign v1 to 2nd col of m1
  REPORT(m1);

  Type m1_sum = m1.sum();              // sum of all elements in m1
  REPORT(m1_sum);

  Type sum_m1 = sum(m1);               // alternative summation
  REPORT(sum_m1);

  Type m1_prod = m1.prod();
  REPORT(m1_prod);

  Type m1_mean = m1.mean();
  REPORT(m1_mean);

  matrix<Type> m2_trans = m2.transpose(); // Transpose of m2
  REPORT(m2_trans);

  vector<Type> m1_diag = m1.diagonal();   // Diagional of m1
  REPORT(m1_diag);

  Type m1_trace = m1.trace();             // Trace of matrix
  REPORT(m1_trace);

  matrix<Type> m3_times_m4 = m3 * m4;   // Matrix product
  REPORT(m3_times_m4);                  // R: m3 %*% m4
                                        
  matrix<Type> m2_times_m2_by_cell = m2.array()*m2.array();
  REPORT(m2_times_m2_by_cell);          // Elementwise product

  matrix<Type> m1_plus_m1_by_cell = m1 + m1; // R: m1 + m1
  REPORT(m1_plus_m1_by_cell);

  vector<Type> m4_colwise_sum = m4.colwise().sum(); // m4 col sums
  REPORT(m4_colwise_sum);

  vector<Type> m4_rowwise = m4.rowwise().sum();     // m4 row sums
  REPORT(m4_rowwise);

  // exp() and log()
  matrix<Type> m2_exp = exp(m2.array());  // exp(m1) does not work
  REPORT(m2_exp);

  // Combining vectors, matrices & arrays
  matrix<Type> a1_plus_m1 = a1.matrix() + m1; // a1+m1.array()
  REPORT(a1_plus_m1);                         //   does not work

  vector<Type> m1_times_v1 = m1*v1;     // matrix-vector product
  REPORT(m1_times_v1);                  // R:  m1 %*% v1

 // Extracting parts of matrices

  matrix<Type> m1_block = m1.block(0,0,1,1);
  REPORT(m1_block); // block starting (0,0) taking 1 row & 1 col


// ARRAYS ============================================================

  int a1_size = a1.size();   // number of elements in array (not dim.)
  REPORT(a1_size);

  vector<int> a1_dim = a1.dim; // dimension of a1
  REPORT(a1_dim);

  int a1_rows = a1.rows();  // Number of rows
  REPORT(a1_rows);

  int a1_cols = a1.cols();  // Number of columns
  REPORT(a1_cols);

  Type a1_11 = a1(0,0);
  REPORT(a1_11);

  a1.col(1) = v1;  // assign v1 to 2nd col of a1
  REPORT(a1);

  Type a1_sum = a1.sum(); // sum of all cells in a1
  REPORT(a1_sum);

  Type sum_a1 = sum(a1); // another way to sum all cells in a1
  REPORT(sum_a1);

  Type a1_mean = a1.mean();           // mean of all cells in a1
  REPORT(a1_mean);

  array<Type> a2_transpose = a2.transpose(); // transpose array a2
  REPORT(a2_transpose);

  array<Type> a2_times_a2_times_a2_by_cell = a2*a2*a2;
  REPORT(a2_times_a2_times_a2_by_cell);

  // TMB loses dim attribute of array when using exp or log;
  //   need to explicitly define dimension attribute of array
  array<Type> a2_exp(a2.dim);
  a2_exp = exp(a2);
  REPORT(a2_exp);

  // Assignment from empty array (used to crash)
  DATA_ARRAY(a0);
  array<Type> a_cpy;
  a_cpy = a0;

  // Splitting a vector.
  // Similar to 'split(c3, f3)' from R
  DATA_FACTOR(f3); // Same length as 'v3'
  vector<vector<Type> > spl3 = split(v3, f3);
  REPORT(spl3);

  return 0;
}
