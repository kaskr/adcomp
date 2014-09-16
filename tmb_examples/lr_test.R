## Illustrate map feature of TMB to perform likelihood ratio tests on a ragged array dataset.
library(TMB)

ngroup <- 5
nrep <- c(5,8,11,13,2)  ## Number of samples per group
mu <- rep(0,ngroup)     ## Mean value per group
sd <- c(1,1,1,2,2)      ## Standard deviation per group

## Simulate data
set.seed(123)
raggedArray <- lapply(1:ngroup,function(i)rnorm(nrep[i],mu[i],sd[i]))

## Prepare data for TMB (ragged array not available):
obs <- unlist(raggedArray)
group <- factor( rep(1:length(raggedArray),sapply(raggedArray,length)) )

compile("lr_test.cpp")
dyn.load(dynlib("lr_test"))

## Both mu's and sd's un-restricted.
full.model <- MakeADFun(data=list(obs=obs,group=group),
                        parameters=list(mu=rep(0,ngroup),sd=rep(1,ngroup)),
                        DLL="lr_test"
                        )

## mu's restricted to be equal
map <- list(mu=factor(c(1,1,1,1,1)))
restricted.model1 <- MakeADFun(data=list(obs=obs,group=group),
                               parameters=list(mu=rep(0,ngroup),sd=rep(1,ngroup)),
                               DLL="lr_test",
                               map=map
                               )

##Both mu's and sd's restricted to be equal
map <- list(mu=factor(c(1,1,1,1,1)),sd=factor(c(1,1,1,1,1)))
restricted.model2 <- MakeADFun(data=list(obs=obs,group=group),
                               parameters=list(mu=rep(0,ngroup),sd=rep(1,ngroup)),
                               DLL="lr_test",
                               map=map
                               )

## Run models:
opt <- do.call("optim",full.model)
opt1 <- do.call("optim",restricted.model1)
opt2 <- do.call("optim",restricted.model2)

LRtest <- function(full,restricted){
    statistic <- 2*(restricted$value - full$value)
    df <- length(full$par) - length(restricted$par)
    p.value <- 1-pchisq(statistic,df=df)
    data.frame(statistic,df,p.value)
}

rbind(LRtest(opt,opt1),
      LRtest(opt1,opt2))
