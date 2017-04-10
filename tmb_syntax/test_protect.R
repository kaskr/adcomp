data <- list()
parameters <- list(p=0)

require(TMB)
compile("test_protect.cpp")
dyn.load(dynlib("test_protect"))

obj <- MakeADFun(data, parameters, DLL="test_protect")
gctorture()
obj$report()
rep <- obj$report(c(p=1))
gctorture(FALSE)
