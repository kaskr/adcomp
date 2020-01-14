#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

if [ ! -f ~/.R/Makevars ]; then
    mkdir -p ~/.R
    echo "CXX = g++ -std=c++11 -Wall -pedantic  " > ~/.R/Makevars
fi

make cran-version
make cran-check
make install
echo "library(TMB);precompile();runExample(all=TRUE)" | R --slave
