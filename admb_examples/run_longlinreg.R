example <- "longlinreg"
example <- sub("\\..*","",example)
print(example)
file <- paste(example,"README.md",sep="/")
if(file.exists(file)) readme <- readLines(file) else readme <- ""
run <- paste(grep(paste0("^\\./",example),readme,value=TRUE),"")
run <- sub("-est","",run)
args <- " -lmn 10"
build <- "cd folder;R --vanilla < longlinreg.R;admb folder.tpl"
run1 <- "cd folder;./folder -args -est"
run2 <- "cd folder;./folder -args"
system(gsub("folder",example,build))
tim1 <- system.time(system(cmd1 <- sub("-args",args,gsub("folder",example,run1))))
tim2 <- system.time(system(cmd2 <- sub("-args",args,gsub("folder",example,run2))))
rep <- try(read.table(paste0(example,"/",example,".std"),header=TRUE))

## Hack: insert par file estimate in std
par <- try(readLines(paste0(example,"/",example,".par")))
par <- grep("^[^#]",par,value=TRUE)
par <- sub("^ ","",par)
par <- do.call("paste",as.list(par))
par <- as.numeric(strsplit(par," ")[[1]])
rep$value[1:length(par)] <- par

save(tim1,tim2,cmd1,cmd2,rep,file=paste0(example,".output.RData"))
