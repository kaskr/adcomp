data <- list(y=1.2)
parameters <- list(mu=1.1, phi=1.1, p=1.1)

require(TMB)
library(tweedie)
compile('tweedie.cpp')
dyn.load(dynlib('tweedie'))

################################################################################

set.seed(1001)
data$y <- rtweedie(1000, mu=2, phi=2, p=1.5)
model <- MakeADFun(data, parameters)
model$fn()
model$gr()
system.time( opt <- nlminb(model$par, model$fn, model$gr) )
sdreport(model)
