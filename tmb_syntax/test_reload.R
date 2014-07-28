## This script tests if re-loading works as intended; i.e. if
## changes apply when editing, re-compiling and reloading a model.
## - If it works the script should display "TRUE".
## - The result currently depends on compiler and operating system.
## - Compiler can be changed to e.g. clang by adding a line "CXX=clang++"
##   to a file named '~/.R/Makevars'

library(TMB)
makeTemplate <- function(output=0){
    txt <-
        "#include <TMB.hpp>
     template<class Type>
     Type objective_function<Type>::operator() (){PARAMETER(a);return output;}"
    txt <- sub("output",output,txt)
    writeLines(txt,"mymodel.cpp")
}

## Make function that returns "1" and compile.
cat("\n")
cat("==============================\n")
cat("Make function that returns '1'\n")
cat("==============================\n")
makeTemplate(1)
compile("mymodel.cpp")
dyn.load(dynlib("mymodel"))
obj <- MakeADFun(list(),list(a=0))
obj$fn() == 1  ## OK ?

## Modify function to return "2" and check output.
## Note: It should not be necessary to explicitly call
## dyn.unload as this is automatically called when
## re-loading a library.
## Library unloading causes TMB to cleanup all external
## pointers before actually unloading. If cleanup fails,
## the library is not unloaded.
cat("\n")
cat("==============================================\n")
cat("Modify function to return '2' and check output\n")
cat("==============================================\n")
makeTemplate(2)
compile("mymodel.cpp")
dyn.load(dynlib("mymodel"))
obj <- MakeADFun(list(),list(a=0))
ok <- obj$fn() == 2  ## OK ?

cat("\n")
cat("OK?",ok,"\n")
