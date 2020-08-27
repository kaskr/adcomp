## ===============================================================
## Profile using intels vtune application
## Script input
## example: Name of example.
example <- Sys.getenv("example")
type <- Sys.getenv("type")
## ===============================================================

Rexe <- paste(Sys.getenv("R_HOME"),"bin/exec/R",sep="/")

## amplxe-cl -help collect
amplxe.options <- c(
    "advanced-hotspots",
    "concurrency",
    "cpugpu-concurrency",
    "disk-io",
    "general-exploration",
    "gpu-hotspots",
    "hotspots",
    "hpc-performance",
    "locksandwaits",
    "memory-access",
    "sgx-hotspots",
    "tsx-exploration",
    "tsx-hotspots")

## inspxe-cl -help collect
##    mi1   Detect Leaks
##    mi2   Detect Memory Problems
##    mi3   Locate Memory Problems
##    ti1   Detect Deadlocks
##    ti2   Detect Deadlocks and Data Races
##    ti3   Locate Deadlocks and Data Races
inspxe.options <- c(
    "mi1",
    "mi2",
    "mi3",
    "ti1",
    "ti2",
    "ti3")

if (!is.na(match(type, amplxe.options))) {
    cmd <- paste("amplxe-cl -collect ",
                 type,
                 " -result-dir ",
                 example,
                 ".profile",
                 " -- ",
                 Rexe,
                 " --vanilla < ",
                 example,".R",sep="")
} else if (!is.na(match(type, inspxe.options))) {
    file.copy(paste0(example, ".R"), paste0(example, ".clean_exit.R") )
    cat("rm(list=ls()); gc()", file=paste0(example, ".clean_exit.R"), append=TRUE)
    cmd <- paste("inspxe-cl -collect=",
                 type,
                 " -result-dir ",
                 example,
                 ".memprofile",
                 " -- ",
                 Rexe,
                 " --vanilla < ",
                 example,".clean_exit.R",sep="")
} else {
    stop("Unknown option collect type: ", type)
}
system(cmd)
