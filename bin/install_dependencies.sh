#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

wget https://cran.r-project.org/src/contrib/Archive/roxygen2/roxygen2_5.0.1.tar.gz
Rscript -e 'pkg <- c("RcppEigen", "numDeriv", "parallel", "knitr", "bookdown", "rsvg", "brew"); if(!all(pkg%in%installed.packages()))install.packages(pkg)'
Rscript -e 'problemPkg <- "roxygen2_5.0.1.tar.gz"; pkg <- c("roxygen2"); if(!all(pkg%in%installed.packages()))install.packages(problemPkg, repos=NULL)'
