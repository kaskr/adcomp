## Experimental RStudio TMB snippet integration
## Gavin Fay & Andrea Havron
## August 16, 2017

setupRSnippets <- function(snipfile="~/.R/snippets/r.snippets",
                           outfile="~/.R/snippets/r.snippets") {

  # read in the text to be added to the snippets
  snips <- readLines(snipfile)
  # add the necessary tabs at the beginnning of each line, 
  # first identifying the headers for each snippet
  headers <- grep("snippet",snips)
  snips[-headers] <- paste("\t ",snips[-headers],sep="")
  #snips <- c(snips[1],paste("\t ",snips[-1],sep=""))
  
    mess <- c("You are about to setup rsnippets for TMB.",
            "Changes will be added to the file:",
            "",
            outfile,
            "")
  invisible(lapply(mess, cat, "\n"))
  ans <- readline("OK? (yes/no) ")
  if(ans == "yes") {
    chk <- grepl(
      snips[headers],
      readLines(outfile)
    )
    if(any(chk))
      message("Skipping because changes seem to have been made already.")
    else {
      #cat(snips, file=outfile, append=TRUE)
      write.table(snips,file=outfile,quote=FALSE,
                  row.names=FALSE,col.names=FALSE,append=TRUE)
      message("Please re-start RStudio for the changes to take place.")
    }
  } else {
    message("Dropping out")
  }
  invisible(NULL)
}

#setupRSnippets(snipfile = "~/sandbox/snips.txt")
