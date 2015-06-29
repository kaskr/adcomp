require(TMB)
library(TMB)
compile("multivariate_distributions.cpp")
dyn.load(dynlib("multivariate_distributions"))
obj <- MakeADFun(data=list(),
                 parameters=list(dummy_par = 0)
		 )

