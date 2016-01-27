## Log cpu and memory usage of test example
logpid <- function(){
    script <-
"logpid() {
  while `kill -s 0 $1`; do
    env -i ps -p $1 -o etime= -o pcpu= -o pmem=;
    sleep .1;
  done;
};
echo elapsed pcpu pmem > FILE;
logpid PID >> FILE &\n"
    script <- gsub("PID", Sys.getpid(), script)
    script <- gsub("FILE", paste0(Sys.getenv("example"), ".logpid") , script)
    cat(script)
    system(script)
}
logpid()
example <- Sys.getenv("example")
source(paste0(example,".R"))
