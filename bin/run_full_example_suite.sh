#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

export OMP_NUM_THREADS=2
mkdir ~/.R
echo "CXX = g++ -Wall -pedantic -Werror" > ~/.R/Makevars
make test
