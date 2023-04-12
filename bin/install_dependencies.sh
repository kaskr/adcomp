#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

Rscript -e 'pkg <- c("RcppEigen", "numDeriv", "parallel", "knitr", "bookdown", "rsvg", "brew", "roxygen2"); if(!all(pkg%in%installed.packages()))install.packages(pkg)'
