#Noncentral chi square MLE
# Compile c++ code and load into R
library(TMB)
compile("dbessel_I_nchisq.cpp")
dyn.load(dynlib("dbessel_I_nchisq"))

n <- 1000
df = 3
ncp = 3.5
data<-rchisq(n, df, ncp)

obj<-MakeADFun(data=list(x=data),parameters=list(k=1,lambda=1))
opt<-nlminb(obj$par,obj$fn,obj$gr,obj$he)
rep<-sdreport(obj)
summary(rep)