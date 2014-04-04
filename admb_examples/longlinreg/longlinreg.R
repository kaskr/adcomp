options("scipen"=100, "digits"=12)
set.seed(123)
nobs<-1000000
x<-seq(0,10, length=nobs)
cat(nobs,"\n",file="longlinreg.dat")
write.table(cbind(2*x+1+rnorm(nobs),x), file="longlinreg.dat", row.names=FALSE, col.names=FALSE, append=TRUE)
