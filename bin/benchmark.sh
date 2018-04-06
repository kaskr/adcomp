#!/bin/bash

THISDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CURRENT_BRANCH=`git rev-parse --abbrev-ref HEAD`
TMB_EXAMPLES=${THISDIR}/../tmb_examples
BENCHMARK_DIR=${THISDIR}/../benchmark/${CURRENT_BRANCH}

mkdir -p ${BENCHMARK_DIR}
rm -rf ${BENCHMARK_DIR}/*
cd ${TMB_EXAMPLES}
git clean -xdf
cp -r ${TMB_EXAMPLES}/* ${BENCHMARK_DIR}

cd ${THISDIR}/..
make install-metis-full
echo "TMB:::precompile()" | R --slave

cd ${BENCHMARK_DIR}
make clean
make test
make logpid_all
