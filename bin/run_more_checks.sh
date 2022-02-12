#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

export OMP_NUM_THREADS=2

if [ ! -f ~/.R/Makevars ]; then
    mkdir -p ~/.R
    echo "CXX = g++ -std=c++11 -Wall -pedantic  " > ~/.R/Makevars
fi

make test-tmb_syntax

