require(TMB)
library(TMB)
compile("matrix_arrays.cpp",CXX="clang++")
dyn.load("matrix_arrays.so")
obj <- MakeADFun(data=list(i=1),
                 parameters=list(p = 0)
		 )

