## Copyright (C) 2013-2015 Kasper Kristensen
## License: GPL-2

##' Compile and run a test example (\code{runExample()} shows all available examples).
##'
##' @title Run one of the test examples.
##' @param name Character name of example.
##' @param all Run all the test examples?
##' @param thisR Run inside this R?
##' @param clean Cleanup before compile?
##' @param exfolder Alternative folder with examples.
##' @param dontrun Build only (don't run) and remove temporary object files ?
##' @param subarch Build in sub-architecture specific folder ?
##' @param ... Passed to \code{\link{compile}}.
runExample <- function(name=NULL,all=FALSE,thisR=TRUE,
                       clean=FALSE,exfolder=NULL,
                       dontrun=FALSE,
                       subarch=TRUE,...){
  cwd <- getwd()
  on.exit(setwd(cwd))
  if(is.null(exfolder))exfolder <- system.file("examples",package="TMB")
  setwd(exfolder)
  arch <- Sys.getenv("R_ARCH")
  if(arch != "" && subarch){
    arch <- sub("/", "", arch)
    if( !file.exists(arch) ){
      dir.create(arch)
      file.copy(dir(pattern="*.[cpp|R]"), arch)
    }
    setwd(arch)
  }
  validExamples <- function(){
    f1 <- sub("\\.[^\\.]*$","",dir(pattern=".R$"))
    f2 <- sub("\\.[^\\.]*$","",dir(pattern=".cpp$"))
    intersect(f1,f2)
  }
  exnames <- validExamples()
  cppnames <- paste0(exnames,".cpp")
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
    lapply(exnames,runExample,thisR=thisR,clean=clean,...)
    return(NULL)
  }
  if(is.null(name)){
    txt <- paste("Examples in " ,"\'",exfolder,"\':","\n\n",sep="")
    cat(txt)
    writeLines(info)
    return(invisible(exnames))
  }
  if(clean){
    cat("Cleanup:\n")
    file.remove(dynlib(name))
    file.remove(paste0(name,".o"))
  }
  if(!file.exists(dynlib(name))){
    cat("Building example",name,"\n")
    time <- system.time(compile(paste0(name,".cpp"),...))
    cat("Build time",time["elapsed"],"seconds\n\n")
    if(dontrun)file.remove(paste0(name,".o"))
  }
  if(!dontrun){
    cat("Running example",name,"\n")
    if(!thisR){
      system(paste("R --vanilla < ",name,".R",sep=""))
    } else {
      source(paste(name,".R",sep=""),echo=TRUE)
    }
  }
}
