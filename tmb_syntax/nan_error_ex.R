# Illustrates how to make the debugger catch a floating point error.
library(TMB)
compile("nan_error_ex.cpp")
dyn.load(dynlib("nan_error_ex"))

## Fit model
obj <- MakeADFun(data=list(lambda=1.5),parameters=list(x=1))
