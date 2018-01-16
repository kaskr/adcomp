#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

if [ ! -f ~/.R/Makevars ]; then
    mkdir -p ~/.R
    echo "CXX = g++ -Wall -pedantic -Werror" > ~/.R/Makevars
fi

make cran-version
make cran-check
make install
echo "library(TMB);precompile();runExample(all=TRUE)" | R --slave
