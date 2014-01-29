## ===============================================================
## Profile using intels vtune application
## Script input
## example: Name of example.
example <- Sys.getenv("example")
## ===============================================================

Rexe <- paste(Sys.getenv("R_HOME"),"bin/exec/R",sep="/")
cmd <- paste("amplxe-cl -collect hotspots -result-dir ",example,".profile"," -- ",Rexe," --vanilla < ",example,".R",sep="")
system(cmd)
cmd <- paste("amplxe-cl -report hotspots -result-dir ",example,".profile"," > ",example,".profile.txt",sep="")
system(cmd)
df <- read.table(paste0(example,".profile.txt"),header=TRUE,sep="\t")
as.matrix(xtabs(CPU.Time~Module,data=df))
