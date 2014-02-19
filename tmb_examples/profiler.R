## ===============================================================
## Profile using intels vtune application
## Script input
## example: Name of example.
example <- Sys.getenv("example")
## ===============================================================

Rexe <- paste(Sys.getenv("R_HOME"),"bin/exec/R",sep="/")

cmd <- paste("amplxe-cl -collect hotspots -result-dir ",example,".profile"," -- ",Rexe," --vanilla < ",example,".R",sep="")
system(cmd)

## Read total duration of profile
cmd <- paste("amplxe-cl -report summary -result-dir ",example,".profile"," | grep Elapsed",sep="")
out <- system(cmd,TRUE,FALSE,TRUE)
endTime <- as.numeric(gsub("[^0-9.]*","",out))
## Bisect to find relevant profile time interval
f <- function(a,b){
  time.filter <- paste(a,b,sep=":")
  cmd <- paste("amplxe-cl -report top-down -result-dir ",example,".profile"," -time-filter ",time.filter,sep="")
  system(cmd,TRUE,FALSE,TRUE)
}
g <- function(x,regexp,tol=0.05,direction=1){
  empty <- length(grep(regexp,f(x[1],x[2])))==0
  dx <- x[2]-x[1]
  print(x)
  if(dx<tol)return(x)
  if(empty){
    x <- x+dx*direction
  } else {
    if(direction>0)
      x <- c(min(x),min(x)+dx/2)
    else
      x <- c(min(x)+dx/2,max(x))
  }
  return(g(x,regexp,tol,direction))
}
if(length(grep("CppAD",f(0,endTime)))!=0){
  begin <- min(g(c(0,2),"CppAD"))
  end <- max(g(c(endTime-1,endTime),"CppAD",direction=-1))
  time.filter <- paste(begin,end,sep=":")
  cat("Profiling time interval",time.filter,"\n")
  cmd <- paste("amplxe-cl -report top-down -result-dir ",example,".profile"," -time-filter ",time.filter," > ",example,".profile.txt",sep="")
  system(cmd)
}
