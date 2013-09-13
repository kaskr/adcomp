##' Compile and run a test example (\code{runExample()} shows all available examples).
##'
##' @title Run one of the test examples.
##' @param name Character name of example.
##' @param all Run all the test examples?
##' @param thisR Run inside this R?
runExample <- function(name=NULL,all=FALSE,thisR=FALSE){
  cwd <- getwd()
  on.exit(setwd(cwd))
  exfolder <- system.file("examples",package="RcppAD")
  setwd(exfolder)
  cppnames <- dir(exfolder,pattern=".cpp$")
  exnames <- sub(".cpp","",cppnames)
  ## Format info as text
  M <- max(nchar(exnames))
  info <- sapply(cppnames,function(x){readLines(x)[[1]]})
  info[substring(info,1,2)!="//"] <- ""
  info <- sub("^// *","",info)
  tmp <- gsub(" ","@",format(paste("\"",exnames,"\"",":",sep=""),width=M+4))
  info <- paste(tmp,info,sep="")
  info <- strwrap(info, width = 60, exdent = M+4)
  info <- gsub("@"," ",info)
  if(all){
    lapply(exnames,runExample,thisR=thisR)
    return(NULL)
  }
  if(is.null(name)){
    txt <- paste("Examples in " ,"\'",exfolder,"\':","\n\n",sep="")
    cat(txt)
    writeLines(info)
    return(invisible(exnames))
  }
  cat("Building example",name,"\n")
  time <- system.time(compile(paste(name,".cpp",sep="")))
  cat("Build time",time["elapsed"],"seconds\n\n")
  cat("Running example",name,"\n")
  if(!thisR){
    system(paste("R --vanilla < ",name,".R",sep=""))
  } else {
    source(paste(name,".R",sep=""),echo=TRUE)
  }
}
