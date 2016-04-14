# Illustrate use of vector, matrix and array operations

library(TMB)
compile("matrix_arrays.cpp")
dyn.load(dynlib("matrix_arrays"))

v1 = c(3,5)
m1 = matrix(c(3,9,2,4),nrow=2)
a1 = matrix(c(1,5,6,7),nrow=2)
    
obj <- MakeADFun(data=list(v1=v1,m1=m1,a1=a1),parameters=list(p = 0),DLL="matrix_arrays")

# Print objects (in alphabetic order) that have been returned by REPORT()
print(obj$env$report())

