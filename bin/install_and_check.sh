#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

mkdir ~/.R
echo "CXX = g++ -Wall -pedantic -Werror" > ~/.R/Makevars
make cran-version
make cran-check
make install
echo "library(TMB);precompile();runExample(all=TRUE)" | R --slave
