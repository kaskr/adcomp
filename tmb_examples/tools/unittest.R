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
  compile <- addHook(TMB::compile,timer=TRUE)
  oneStepPredict <- addHook(TMB::oneStepPredict,timer=TRUE,result=TRUE)
  
  setwd(dirname(example))
  source(paste0(basename(example), ".R"), echo=TRUE)

  ## Strip off large irrelevant output
  strip.off <- function(x) {
      if(is.list(x)) {
          x[] <- lapply(unclass(x), strip.off)
      }
      else if (is.environment(x)) {
          x <- new.env() ## empty
      } else {
          ## Do nothing
      }
      x
  }
  .results <- strip.off(.results)

  if(!file.exists(paste0(basename(example),".expected.RData"))){
    outfile <- paste0(basename(example),".expected.RData")
    save(.timings,.results,file=outfile)
  }
  outfile <- paste0(basename(example),".output.RData")
  save(.timings,.results,file=outfile)
  
} else {
  report.level <- Sys.getenv("report_level")
  report.full <- (report.level == "1")
  ## Report of diffs
  f1 <- dir(pattern = ".expected.RData$", recursive=TRUE)
  f2 <- sub("\\.expected\\.RData$","\\.output\\.RData",f1)
  report <- function(f1,f2,full.timings=FALSE,full.diff=FALSE){
    if(!(file.exists(f1)&file.exists(f2)))return(c("NA"=NA))
    diff <- function(x,y){
      if(is.list(x)&is.list(y)){
          x <- unclass(x) ## avoid as.list.sdreport
          y <- unclass(y) ## avoid as.list.sdreport
          ## Allow to compare with old expected output (that did not
          ## have 'env' as part of sdreport output):
          keep <- function(x)
              !is.environment(x)
          keepx <- sapply(x,keep)
          keepy <- sapply(y,keep)
          Map(diff, x[keepx], y[keepy])
      }
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
    ## Total runtime (exclude TMB::compile)
    totaltime <- function(x){
        x <- x[names(x) != "TMB::compile"]
        sum(sapply(x,function(x)x["elapsed"]))
    }
    t1 <- totaltime(e1$.timings)
    t2 <- totaltime(e2$.timings)
    ## Total compilation time
    compilationtime <- function(x){
        x <- x[names(x) == "TMB::compile"]
        sum(sapply(x,function(x)x["elapsed"]))
    }
    tc1 <- compilationtime(e1$.timings)
    tc2 <- compilationtime(e2$.timings)
    c(summary(as.numeric(d))[c("Min.","Median","Max.")],  ## Results
      totaltime=t2, timeindex=t2/t1,                      ## Runtime
      ctime=tc2, ctimeindex=tc2/tc1)                      ## Compilation time
  }
  sink("REPORT.md")
  options(width = 100)
  cat("Example overview:\n-----------------\n")
  runExample(exfolder=".",subarch=FALSE)
  cat("\n")
  mat <- do.call("rbind",Map(report,f1,f2))
  rownames(mat) <- sub(".expected.RData","",rownames(mat))
  if (TRUE) {
      cat("\nResults (absolute error):\n-------------------------\n")
      print(mat[ , c("Min.", "Median", "Max.")])
  }
  if (TRUE) {
      cat("\nTimings:\n--------\n")
      print(mat[ , c("totaltime", "timeindex")])
  }
  if (report.full) {
      cat("\nCompilation time:\n-----------------\n")
      print(mat[ , c("ctime", "ctimeindex")])
  }
  if (report.full) cat("\nTiming details:\n---------------\n")
  res <- Map(report,f1,f2,full.timings=TRUE)
  Example <- sub(".expected.RData","",rep(names(res),sapply(res,length)))
  Function <- sub(".*::(.*).elapsed","\\1",unlist(lapply(res,names)))
  tab <- xtabs(unlist(res)~Example+Function)
  names(dimnames(tab)) <- NULL
  if (report.full) print(tab)
  if (report.full) cat("\nResult details:\n---------------\n")
  res <- Map(report,f1,f2,full.diff=TRUE)
  Example <- sub(".expected.RData","",rep(names(res),sapply(res,length)))
  Function <- sub(".*::(.*)","\\1",unlist(lapply(res,names)))
  Function <- formatC(Function,width=max(nchar(Function)))
  tab <- xtabs(unlist(res)~Example+Function)
  names(dimnames(tab)) <- NULL
  if (report.full) print(tab)
  ## Tolerances (default = 1e-8 - modify for selected slots)
  abs.tol <- structure(rep(1e-8, ncol(tab)),
                       .Names=gsub(" ","",colnames(tab)))
  abs.tol["nlminb.objective"] <- 1e-04
  abs.tol["nlminb.par"] <- 1e-05
  abs.tol["optim.par"] <- 1e-06
  abs.tol["optim.value"] <- 1e-07
  abs.tol["sdreport.cov.fixed"] <- 1e-04
  abs.tol["sdreport.diag.cov.random"] <- 1e-07
  abs.tol["sdreport.gradient.fixed"] <- 6e-04
  abs.tol["sdreport.par.fixed"] <- 1e-05
  abs.tol["sdreport.par.random"] <- 1e-06
  abs.tol["sdreport.value"] <- 1e-06
  abs.tol <- abs.tol[gsub(" ","",colnames(tab))]
  passed <- t(t(tab) < abs.tol)
  all.passed <- t(t(apply(passed, 1, all)))
  colnames(all.passed) <- "passed"
  cat("\nAccuracy tests:\n---------------\n")
  print(all.passed)
  ## What not passed ?
  if(!all(all.passed)) {
    names(dimnames(tab)) <- c("Example", "Slot")
    df <- as.data.frame(t(tab), responseName="Error")
    df <- df[c("Example", "Slot", "Error")]
    df$tol <- abs.tol
    cat("\nNot passed:\n-----------\n")
    print(df[!(t(tab) < abs.tol),,drop=FALSE])
  }
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
