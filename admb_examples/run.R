example <- Sys.getenv("example")
example <- sub("\\..*","",example)
print(example)
build <- "cd folder;admb -r folder.tpl"
run1 <- "cd folder;./folder -est"
run2 <- "cd folder;./folder"
system(gsub("folder",example,build))
tim1 <- system.time(system(gsub("folder",example,run1)))
tim2 <- system.time(system(gsub("folder",example,run2)))
rep <- read.table(paste0(example,"/",example,".std"),header=TRUE)
save(tim1,tim2,rep,file=paste0(example,".output.RData"))
