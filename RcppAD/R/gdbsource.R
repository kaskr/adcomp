##' Source R-script through gdb to get backtrace.
##'
##' This function is useful for debugging templates.
##' If a script aborts e.g. due to an out-of-bound index operation
##' it should be fast to locate the line that caused the problem by
##' running \code{gdbsource(file)}.
##' @title Source R-script through gdb to get backtrace.
##' @param file Your R script
##' @return Object of class \code{backtrace}
gdbsource <- function(file){
  gdbscript <- tempfile()
  cat("run\nbt\nquit\n",file=gdbscript)
  cmd <- paste("R --vanilla < ",file," -d gdb --debugger-args=\"-x",
               gdbscript,"\"")
  txt <- system(cmd,intern=TRUE,ignore.stdout=FALSE,ignore.stderr=TRUE)
  attr(txt,"file") <- file
  class(txt) <- "backtrace"
  txt
}
print.backtrace <- function(x,...){
  ## Backtrace begins here
  line <- grep("#0",x)
  ## Assume cpp file same name as r file
  pattern <- gsub("[R|r]$","",attr(x,"file"))
  x <- x[line:length(x)]
  x <- grep(pattern,x,value=TRUE)
  cat(paste(x,"\n"))  
}
