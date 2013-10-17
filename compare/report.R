require(RcppAD)
f1 <- dir("../rcppad_examples",".expected.RData")
f2 <- dir("../admb_examples",".output.RData")
examples <- intersect(sub("\\..*","",f1),sub("\\..*","",f2))
compare <- function(example){
  owd <- getwd()
  on.exit(setwd(owd))
  setwd("../admb_examples")
  load(paste0(example,".output.RData"))
  setwd("../rcppad_examples")
  load(paste0(example,".expected.RData"))
  ok <- all(as.character(rep$name)==rownames(sum))
  if(!ok){
    cat("Example:",example,"\n")
    cat("Sdreport names must match between the two examples.\n")
    stop()
  }
  sum <- summary(.results$`RcppAD::sdreport`)
  diff <- sum-rep[,3:4]
  res <- apply(abs(diff),2,max)
  names(res) <- c("Max norm est-diff","Max norm sd-diff")
  tim <- as.numeric(.timings$`RcppAD::sdreport`["elapsed"])
  tot <- sum(sapply(.timings,function(x)x["elapsed"]))
  res2 <- c("rcppad est time"=tot-tim,"rcppad sdrep time"=tim)
  res3 <- c("admb est time"=as.numeric(tim1["elapsed"]),
            "admb sdrep time"=as.numeric(tim2["elapsed"]-tim1["elapsed"]))
  c(res,res2,res3)
}
t(sapply(examples,compare))

