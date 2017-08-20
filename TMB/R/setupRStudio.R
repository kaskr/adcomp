## Experimental RStudio integration
setupRStudio <- function(file = "~/.Rprofile",
                         snipRfile = "~/.R/snippets/r.snippets",
                         snipCppfile = "~/.R/snippets/c_cpp.snippets") {
    on.exit(
        message("Please re-start RStudio for the changes to take place.")
    )
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
                ## Auto completion needs TMB and Eigen on system includes
                if (.Platform$OS.type=="windows") {
                  ## Overload system.file
                  system.file <- function(...){
                    ans <- base::system.file(...)
                    chartr("\\\\", "/", shortPathName(ans))
                  }
                }
                definc <- Sys.getenv("CPLUS_INCLUDE_PATH")
                tmbinc <- system.file("include", package="TMB")
                eiginc <- system.file("include", package="RcppEigen")
                inc <- c(definc, tmbinc, eiginc)
                inc <- paste(inc[inc != ""], collapse=.Platform$path.sep)
                Sys.setenv(CPLUS_INCLUDE_PATH = inc)
            }
        } )
'
    mess <- c("You are about to setup Rstudio with TMB.",
              "Changes will be added to the file:",
              "",
              file,
              "")
    invisible(lapply(mess, cat, "\n"))
    getYesOrNo <- function() {
        repeat {
            ans <- readline("OK? (yes/no) ")
            if (ans %in% c("yes", "no")) break;
            message("Please say 'yes' or 'no'")
        }
        ans
    }
    ans <- getYesOrNo()
    if(ans == "yes") {
        ## Create ~/.Rprofile if not exists
        if (!file.exists(file)) file.create(file)
        ## Read ~/.Rprofile and remove change if previously made
        oldlines <- readLines(file)
        codelines <- strsplit(code,"\n")[[1]][-1]
        begin <- which(head(codelines, 1) == oldlines)[1]
        if (!is.na(begin)) {
            if (begin > 1) begin <- begin - 1
            end <- which(tail(codelines, 1) == oldlines)
            end <- min(end[end>begin])
            oldlines <- oldlines[-(begin:end)]
            message("Removing old changes from ", file)
            writeLines(oldlines, file)
        }
        message("Adding changes to ", file)
        cat(code, file=file, append=TRUE)
    }

    ## Experimental RStudio TMB snippet integration
    ## Gavin Fay & Andrea Havron
    ##
    rsnips <- getRsnips()
    headers <- grep("snippet",rsnips)
    rsnips[-headers] <- paste0("\t",rsnips[-headers])
    cppsnips <- getCppsnips()
    cheaders <- grep("snippet",cppsnips)
    cppsnips[-cheaders] <- paste0("\t",cppsnips[-cheaders])

    mess <- c("",
              "You are about to setup snippets for TMB.",
              "Changes will be added to the files:",
              "",
              snipRfile,
              snipCppfile,
              "")
    invisible(lapply(mess, cat, "\n"))
    ans <- getYesOrNo()
    if(ans == "yes") {
        if (file.exists(snipRfile) &&
            any( grepl(rsnips[headers][1], readLines(snipRfile)) ) ) {
            message("Skipping because changes seem to have been made already.")
        }
        else {
            dir.create(dirname(snipRfile), showWarnings=FALSE, recursive=TRUE)
            if (!file.exists(snipRfile)) file.create(snipRfile)
            cat(paste(rsnips, collapse="\n"), file=snipRfile, append=TRUE)
            dir.create(dirname(snipCppfile), showWarnings=FALSE, recursive=TRUE)
            if(!file.exists(snipCppfile))file.create(snipCppfile)
            cat(paste(cppsnips, collapse="\n"), file=snipCppfile, append=TRUE)
        }
    }
    invisible(NULL)
}

## R snippets, paste as text, replace \$ with \\$
getRsnips <- function() {
    snips <-
'
snippet tmb.template
## Load TMB `r require(TMB)`
library(TMB)

## Make C++ file
TMB::template("${1:model_name}.cpp")

## Compile and load the model
compile("${1:model_name}.cpp")
dyn.load(dynlib("${1:model_name}"))

## Data and parameters
data <- list(x=rivers)
parameters <- list(mu=0, logSigma=0)

## Make a function object
obj <- MakeADFun(data, parameters, DLL="${1:model_name}")

## Call function minimizer
opt <- nlminb(obj\\$par, obj\\$fn, obj\\$gr)

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr

'
    strsplit(snips, "\n")[[1]][-1]
}


## R snippets, paste as text, replace \$ with \\$
## FIXME: C++ snippets doesn't seem to work
getCppsnips <- function() {
    snips <-
'
snippet tmb.template
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);
  PARAMETER(mu);
  PARAMETER(logSigma);

  Type f = 0;
  f -= dnorm(x, mu, exp(logSigma), true).sum();

  return f;
}

'
    strsplit(snips, "\n")[[1]][-1]
}
