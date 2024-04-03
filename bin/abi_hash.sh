#!/bin/bash -e
TMPFILE="$( mktemp --suffix=.cpp )"
git show $1:TMB/inst/include/tmb_core.hpp > ${TMPFILE}
echo 'TMB_CALLDEFS'            >> ${TMPFILE}
echo 'TMB_CALLABLES(PKG)'      >> ${TMPFILE}
ABI="$( gcc -E ${TMPFILE} | tail -n 2 )"
rm ${TMPFILE}
echo ${ABI} | md5sum
