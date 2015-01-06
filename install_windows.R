if(.Platform$OS.type != "windows")stop("This script is only for windows")

## ===== Default locations - change for non-standard setup
fileLocations <- function(){
  rhome <- Sys.getenv("R_HOME")
  rbin <- paste0(rhome,"/bin/",.Platform$r_arch)
  rtools <- c("C:/Rtools",
              "~/Rtools",
              paste0(rhome,"/../../Rtools"),
              paste0(rhome,"/../Rtools") )
  gccbin <- paste0(rtools,"/gcc-4.6.3/bin")
  gdbbin <- paste0(rtools,"/gcc-4.6.3/bin64")
  ans <- list(rhome=rhome,rbin=rbin,rtools=rtools,gccbin=gccbin,gdbbin=gdbbin)
  lapply(ans,function(x)x[file.exists(x)][1])
}

installRtools <-function(){
  owd <- getwd()
  fl <- fileLocations()
  if(!is.na(fl$rtools)){
    cat("Rtools already installed\n")
    return(NULL)
  }
  ## ===== Download and run Rtools installer
  rtools <- "http://cran.r-project.org/bin/windows/Rtools/Rtools31.exe"
  dest <- paste(tempdir(),"Rtools31.exe",sep="\\")
  download.file(rtools,dest,mode="wb")
  setwd(tempdir())
  shell("Rtools31.exe")
  setwd(owd)
}

setPath <- function(...){
  fl <- fileLocations()
  ## ===== Add R to PATH  
  path <- Sys.getenv("PATH")
  path <- paste(path.expand(fl$rbin),path,sep=";")
  ## ===== Add gcc to PATH
  if(is.na(fl$rtools))stop("Failed to locate Rtools folder")
  path <- paste(path.expand(fl$gccbin),path,sep=";")
  ## ===== Add gdb to PATH
  path <- paste(path.expand(fl$gdbbin),path,sep=";")
  ## ===== Add unix tools to PATH
  cwpath <- paste0(path.expand(fl$rtools),"/bin")
  path <- paste(cwpath,path,sep=";")
  ## ====== Update PATH
  Sys.setenv(PATH=path)
  ## ====== Disable cygwin warnings
  Sys.setenv(CYGWIN="nodosfilewarning")
}

createMakevars <- function(){
  file <- "~/.R/Makevars"
  if(file.exists(file))return(NULL)
  dir.create("~/.R")
  rhome <- Sys.getenv("R_HOME")
  p <- paste0(rhome,"/etc/x64/Makeconf")
  cat("Changing Makevars from\n")
  x <- grep("-Wall",readLines(p),val=TRUE)
  print(x)
  x <- sub("-Wall ","",x)
  x <- c(x,"PKG_CPPFLAGS += -DSUPPORT_OPENMP")
  cat("to\n")
  print(x)
  writeLines(x,file)
}

modifyPathOnAttach <- function(){
  winpath <- c("fileLocations<-",deparse(fileLocations),
               "setPath<-",deparse(setPath),
               ".onAttach<-function(libname, pkgname){setPath()}")
  writeLines(winpath,"TMB/R/winpath.R")
}

if(!("TMB" %in% installed.packages())){
  installRtools()
  setPath()
  createMakevars()
  modifyPathOnAttach()
  shell("make install")
} else {
  cat("TMB already installed. To reinstall run:\n")
  cat("  remove.packages(\"TMB\")\n")
  cat("  source(\"install_windows.R\")\n")
}
