## Copyright (C) 2013-2015 Kasper Kristensen
## License: GPL-2

##' Source R-script through gdb to get backtrace.
##'
##' This function is useful for debugging templates.
##' If a script aborts e.g. due to an out-of-bound index operation
##' it should be fast to locate the line that caused the problem by
##' running \code{gdbsource(file)}.
##' Alternatively, If more detailed debugging is required,  then
##' \code{gdbsource(file,TRUE)} will provide the full backtrace followed
##' by an interactive gdb session where the individual frames can be inspected.
##' Note that templates should be compiled without optimization and with debug
##' information in order to provide correct line numbers:
##' \itemize{
##' \item On Linux/OS X use \code{compile(cppfile,"-O0 -g")}.
##' \item On Windows use \code{compile(cppfile,"-O1 -g",DLLFLAGS="")} (lower
##' optimization level will cause errors).
##' }
##' @title Source R-script through gdb to get backtrace.
##' @param file Your R script
##' @param interactive Run interactive gdb session?
##' @return Object of class \code{backtrace}
gdbsource <- function(file,interactive=FALSE){
  if(!file.exists(file))stop("File '",file,"' not found")
  if(.Platform$OS.type=="windows"){
    return(.gdbsource.win(file,interactive))
  }
  gdbscript <- tempfile()
  if(interactive){
    gdbcmd <- c(paste("run --vanilla <",file),
                "bt")
    gdbcmd <- paste(gdbcmd,"\n",collapse="")
    cat(gdbcmd,file=gdbscript)
    cmd <- paste("R -d gdb --debugger-args=\"-x",gdbscript,"\"")
    system(cmd,intern=FALSE,ignore.stdout=FALSE,ignore.stderr=TRUE)
    return(NULL)
  } else {
    cat("run\nbt\nquit\n",file=gdbscript)
    cmd <- paste("R --vanilla < ",file," -d gdb --debugger-args=\"-x",
                 gdbscript,"\"")
    txt <- system(cmd,intern=TRUE,ignore.stdout=FALSE,ignore.stderr=TRUE)
    attr(txt,"file") <- file
    class(txt) <- "backtrace"
    return(txt)
  }
}
## Windows case
.gdbsource.win <- function(file,interactive=FALSE){
  gdbscript <- tempfile()
  txt <- paste("set breakpoint pending on\nb abort\nrun --vanilla -f",
               file, "\nbt\n")
  cat(txt, file=gdbscript)
  cmd <- paste("gdb Rterm -x", gdbscript)
  if(interactive){
    cmd <- paste("start",cmd)
    shell(cmd)
    return(NULL)
  }
  else {
    txt <- system(cmd,intern=TRUE,ignore.stdout=FALSE,ignore.stderr=TRUE)
    attr(txt,"file") <- file
    class(txt) <- "backtrace"
    return(txt)
  }
}

##' If \code{gdbsource} is run non-interactively (the default) only
##' the relevant information will be printed. Note that this will only
##' work if the cpp file and the R file share the same base name.
##'
##' @title Print problematic cpp line number.
##' @param x Backtrace from \code{gdbsource}
##' @param ... Not used
##' @rdname gdbsource
##' @method print backtrace
##' @S3method print backtrace
##' @return NULL
print.backtrace <- function(x,...){
  ## Backtrace begins here
  line <- grep("#0",x)
  ## Assume cpp file same name as r file
  pattern <- gsub("[R|r]$","",attr(x,"file"))
  x <- x[line:length(x)]
  x <- grep(pattern,x,value=TRUE)
  cat(paste(x,"\n"))  
}
