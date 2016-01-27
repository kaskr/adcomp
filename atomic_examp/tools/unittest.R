## ===============================================================
## Script input
## example: Name of example.
## If missing, generate a report of output versus expected output
## for all examples.
example <- Sys.getenv("example")
## ===============================================================
library(TMB)

if(example!=""){
  .timings <- list()
  .results <- list()
  addHook <- function(f,timer=FALSE,result=FALSE){
    name <- deparse(substitute(f))
    g <- function(...){
      tim <- system.time(ans <- f(...))
      if(timer){
        li <- list(tim);names(li) <- name
        .GlobalEnv$.timings <- c(.GlobalEnv$.timings,li)
      }
      if(result){
        li <- list(ans);names(li) <- name
        .GlobalEnv$.results <- c(.GlobalEnv$.results,li)
      }
      ans
    }
    g
  }
  MakeADFun <- addHook(TMB::MakeADFun,timer=TRUE)
  sdreport <- addHook(TMB::sdreport,timer=TRUE,result=TRUE)
  optim <- addHook(stats::optim,timer=TRUE,result=TRUE)
  nlminb <- addHook(stats::nlminb,timer=TRUE,result=TRUE)
  
  runExample(example,exfolder=".",thisR=TRUE)
  
  if(!file.exists(paste0(example,".expected.RData"))){
    outfile <- paste0(example,".expected.RData")
    save(.timings,.results,file=outfile)
  }
  outfile <- paste0(example,".output.RData")
  save(.timings,.results,file=outfile)
  
} else {
  ## Report of diffs
  f1 <- dir(pattern = ".expected.RData$")
  f2 <- sub("\\.expected\\.","\\.output\\.",f1)
  report <- function(f1,f2,full.timings=FALSE,full.diff=FALSE){
    if(!(file.exists(f1)&file.exists(f2)))return(c("NA"=NA))
    diff <- function(x,y){
      if(is.list(x)&is.list(y))Map(diff,x,y)
      else if((!is.integer(x))&(is.numeric(x)|is.matrix(y))&length(x)>0)max(abs(x-y))
      else NULL
    }
    e1 <- local({load(f1);environment()})
    e2 <- local({load(f2);environment()})
    d <- unlist(diff(e1$.results,e2$.results))
    if(full.diff){
      return( d )
    }
    if(full.timings){
      return( sapply(e2$.timings,function(x)x["elapsed"]) )
    }
    totaltime <- function(x)sum(sapply(x,function(x)x["elapsed"]))
    t1 <- totaltime(e1$.timings)
    t2 <- totaltime(e2$.timings)
    c(summary(as.numeric(d))[c("Min.","Median","Max.")],totaltime=t2,timeindex=t2/t1)
  }
  sink("REPORT.md")
  cat("Example overview:\n-----------------\n")
  runExample(exfolder=".")
  cat("\n")
  mat <- do.call("rbind",Map(report,f1,f2))
  rownames(mat) <- sub(".expected.RData","",rownames(mat))
  cat("\nResults and timings:\n--------------------\n")
  print(mat)
  cat("\nTiming details:\n---------------\n")
  res <- Map(report,f1,f2,full.timings=TRUE)
  Example <- sub(".expected.RData","",rep(names(res),sapply(res,length)))
  Function <- sub(".*::(.*).elapsed","\\1",unlist(lapply(res,names)))
  tab <- xtabs(unlist(res)~Example+Function)
  names(dimnames(tab)) <- NULL
  print(tab)
  cat("\nResult details:\n---------------\n")
  res <- Map(report,f1,f2,full.diff=TRUE)
  Example <- sub(".expected.RData","",rep(names(res),sapply(res,length)))
  Function <- sub(".*::(.*)","\\1",unlist(lapply(res,names)))
  Function <- formatC(Function,width=max(nchar(Function)))
  tab <- xtabs(unlist(res)~Example+Function)
  names(dimnames(tab)) <- NULL
  print(tab)
  sink()
  ## Markdown
  md <- function(file){
    li <- readLines(file)
    i <- grep("^---",li)
    i <- c(i-1,i)
    li[-i] <- paste0("    ",li[-i])
    writeLines(li,file)
  }
  md("REPORT.md")
}
