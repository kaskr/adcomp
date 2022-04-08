## Run transform.R with modified template
tmpfile <- tempfile()
newscript <- gsub("transform","transform_parallel",readLines("transform.R"))
newscript[4] <- "openmp(max=TRUE)"
writeLines(newscript,tmpfile)
source(tmpfile,echo=TRUE)

## Compare parallel result with non-parallel result
file.copy("transform.expected.RData","transform_parallel.expected.RData",TRUE)
