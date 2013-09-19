timeSteps=2000


sd=0.5
sdObs=1;



d=numeric(timeSteps)
obs=d;

##init state
d[1] = rnorm(1);
i=1;
obs[i,] = d[i,] + rnorm(1,0,sdObs)
for(i in 2:timeSteps){
    d[i] = d[i-1] + rnorm(1,0,sd)
    obs[i] = d[i] + rnorm(1,0,sdObs)
}

## make ADMB data file
cat(timeSteps,"\n",file="rw.dat");
##cat(stateDim,"\n",file="rw.dat",append=TRUE);
write.table(obs,file="rw.dat",append=TRUE,row.names=FALSE,col.names=FALSE)


t1=system.time(system("./rw -noinit"))

t2=system.time(system("./rw -noinit -nohess"))

t3=system.time(system("./rw -noinit -nohess -ilmn 5"))

t4=system.time(system("./rw -noinit -nohess -lmn2 5"))


library(RcppAD)
setwd("../../rcppad_examples/")

compile("rw1d.cpp")

dyn.load("rw1d.so")

data <- list(obs=obs)
parameters <- list(
  u=data$obs*0,
  logsd=sd*0,
  logsdObs=sdObs*0
  )
newtonOption(smartsearch=FALSE)
obj <- MakeADFun(data,parameters,random="^u",DLL="rw1d")
obj$hessian <- TRUE

obj$fn()
obj$gr()
system.time(opt <- do.call("optim",obj))
system.time(opt2 <- nlminb(obj$par,obj$fn,obj$gr))

pl <- obj$env$parList() ## <-- List of predicted random effects
plot(pl$u)
lines(obs,col=2)

