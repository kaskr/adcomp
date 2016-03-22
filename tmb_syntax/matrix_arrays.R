library(TMB)
compile("matrix_arrays.cpp")
dyn.load(dynlib("matrix_arrays"))
obj <- MakeADFun(data=list(i=1),
                 parameters=list(p = 0)
		 )
