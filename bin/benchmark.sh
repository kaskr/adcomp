#!/bin/bash -e

if [ "$1" == '' ]
then
    echo 'usage: bin/benchmark.sh test_name'
    echo 'Runs the specified test. If test_name is all, all the tests are run'
    exit 1
fi
TEST_NAME="$1"
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
if [ "$TEST_NAME" == 'all' ]
then
    echo "TMB:::precompile()" | R --slave
fi

cd ${BENCHMARK_DIR}
make clean
if [ "$TEST_NAME" == 'all' ]
then
    make test
    make logpid_all
else
    make $TEST_NAME
    make $TEST_NAME.logpid
fi
