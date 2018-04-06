## Run register_atomic.R with modified template
tmpfile <- tempfile()
newscript <- gsub("register_atomic","register_atomic_parallel",readLines("register_atomic.R"))
writeLines(newscript,tmpfile)
source(tmpfile,echo=TRUE)

## Compare parallel result with non-parallel result
file.copy("register_atomic.expected.RData","register_atomic_parallel.expected.RData",TRUE)
