## Experimental RStudio integration
setupRStudio <- function(file="~/.Rprofile") {
code <- 
'
######## TMB - setup RStudio
setHook(packageEvent("TMB", "onLoad"),
        function(...) {
            if("tools:rstudio" %in% search()) {
                tmb.env <- asNamespace("TMB")
                compile.orig <- tmb.env$compile
                unlockBinding("compile", tmb.env)
                ## Rstudio handle compilation errors:
                rs.env <- as.environment("tools:rstudio")
                tmb.env$compile <- function(file,...) {
                    .Call("rs_sourceCppOnBuild",
                          file, FALSE, FALSE)
                    status <- try( compile.orig(file, ...) )
                    succeeded <- (status == 0)
                    .Call("rs_sourceCppOnBuildComplete",
                          succeeded, "")
                    if(!succeeded) stop("Compilation failed")
                    status
                }
                ## Bind "sourceCpp" button to TMB compile
                rcpp.env <- asNamespace("Rcpp")
                unlockBinding("sourceCpp", rcpp.env)
                rcpp.env$sourceCpp <- tmb.env$compile
            }
        } )
'
mess <- c("You are about to setup Rstudio with TMB.",
          "Changes will be added to the file:",
          "",
          file,
          "")
invisible(lapply(mess, cat, "\n"))
ans <- readline("OK? (yes/no) ")
if(ans == "yes") {
    chk <- grepl(
        "######## TMB - setup RStudio",
        readLines(file)
    )
    if(any(chk))
        message("Skipping because changes seem to have been made already.")
    else {
        cat(code, file=file, append=TRUE)
        message("Please re-start RStudio for the changes to take place.")
    }
} else {
    message("Dropping out")
}

# call to the function that adds R snippets. Need a similar one to add C++ snippets
setupRSnippets(snipfile = file.path(system.file("examples",package="TMB"),
                                    "snips.txt"))
# GF put the snipfile in the examples folder but there is probably a better place for it.

invisible(NULL)
}
