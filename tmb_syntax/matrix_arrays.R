# Demonstrates use of vectors, matrices and arrays. 
# Results of operations reported back to R via REPORT()

model_data = list(v1 = c(9,11), 
                  v2 = c(1,2,3,4), 
                  v3 = c(10,20,30,40,50),
                  m1 = matrix(c(1,2,3,4),        nrow=2, ncol=2), 
                  m2 = matrix(c(3,4,5,6,7),      nrow=3, ncol=5),
                  m3 = matrix(c(8,9,10,11),      nrow=1, ncol=4),
                  m4 = matrix(c(10,20,30,40,50), nrow=4, ncol=5),
                  a1 = array(c(1,2,3,4),         dim = c(2,2)),
                  a2 = array(c(8,9,10,11,12),    dim = c(7,5)))

parameters = list(p=1.0)

require(TMB)
compile("matrix_arrays.cpp")
dyn.load(dynlib("matrix_arrays"))

model = MakeADFun(model_data, parameters,DLL="matrix_arrays")
print(model$report())  # Note: order of variables NOT the same as .cpp file


