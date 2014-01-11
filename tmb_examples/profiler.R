## ===============================================================
## Profile using intels vtune application
## Script input
## example: Name of example.
example <- Sys.getenv("example")
## ===============================================================

Rexe <- paste(Sys.getenv("R_HOME"),"bin/exec/R",sep="/")
cmd <- paste("amplxe-cl -collect hotspots -result-dir ",example,".profile"," -- ",Rexe," --vanilla < ",example,".R",sep="")
system(cmd)
