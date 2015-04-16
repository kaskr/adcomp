readadmb <- function(file="nmix.dat"){
  d <- readLines(file)
  i <- grep("#",d)[-1]
  ## Data names
  nm <- gsub(" ","",gsub("#","",d[i]))
  ind <- Map(":",
             i+1,
             c(i[-1]-1,length(d)))
  spl <- lapply(ind,function(x)d[x])
  names(spl) <- nm
  asnum <- function(x){
    x <- as.numeric(unlist(strsplit(x," ")))
    x[!is.na(x)]
  }
  dat <- lapply(spl,asnum)
  dat
}
