#!/bin/bash

THISDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

###########################################################
R --version

###########################################################
## make check

###########################################################
cd ${THISDIR}/../tmb_syntax
R --slave < check_simulations.R > out.txt
if grep -q FALSE out.txt; then
    echo "Error: At least one test not passed"
    cat out.txt
    exit 1
fi

###########################################################
cd ${THISDIR}/../tmb_examples
make all
cat REPORT.md
if grep -q FALSE REPORT.md; then
    echo "Error: At least one test not passed"
    exit 1
fi
