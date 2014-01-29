require(TMB)
f1 <- dir("../tmb_examples",".output.RData")
f2 <- dir("../admb_examples",".output.RData")
examples <- intersect(sub("\\..*","",f1),sub("\\..*","",f2))
compare <- function(example){
  owd <- getwd()
  on.exit(setwd(owd))
  setwd("../admb_examples")
  load(paste0(example,".output.RData"))
  setwd("../tmb_examples")
  load(paste0(example,".output.RData"))
  sum <- summary(.results$`TMB::sdreport`)
  M <- min(nrow(sum),nrow(rep))
  if(M<1)stop("No output")
  rep <- rep[1:M,] ## In case admb example has sdreport - we don't compare that yet
  sum <- sum[1:M,]
  ok <- all(as.character(rep$name)==rownames(sum))
  if(!ok){
    cat("Example:",example,"\n")
    cat("Sdreport names must match between the two examples.\n")
    cat("TMB example:\n")
    print(table(rownames(sum)))
    cat("ADMB example:\n")
    print(table(rep$name))
    stop()
  }
  nfixed <- nrow(summary(.results$`TMB::sdreport`,"fixed")) ## number of fixed effects
  ##diff <- sum-rep[,3:4]
  i <- 1:nfixed
  diff <- sum[i,]-rep[i,3:4] ## Only compare fixed effects
  res <- apply(abs(diff),2,max)
  names(res) <- c("Max norm est-diff","Max norm sd-diff")
  diff <- sum[-i,]-rep[-i,3:4] ## Only compare random effects
  resrf <- apply(abs(diff),2,max)
  names(resrf) <- c("Max norm rfest-diff","Max norm rfsd-diff")
  tim <- as.numeric(.timings$`TMB::sdreport`["elapsed"])
  tot <- sum(sapply(.timings,function(x)x["elapsed"]))
  res2 <- c("tmb est time"=tot-tim,"tmb sdrep time"=tim)
  res3 <- c("admb est time"=as.numeric(tim1["elapsed"]),
            "admb sdrep time"=as.numeric(tim2["elapsed"]-tim1["elapsed"]))
  c(res,resrf,res2,res3)
}
t(sapply(examples,compare))

