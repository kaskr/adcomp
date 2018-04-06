#!/bin/bash -e
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
    echo $*
    eval $*
}
# -----------------------------------------------------------------------------
THISDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CPPAD_DIR=${THISDIR}/../../CppAD
TMB_INCLUDE=${THISDIR}/../TMB/inst/include
BUILD_CPPAD=${CPPAD_DIR}/build_cppad_tmb
# -----------------------------------------------------------------------------
if [ ! -d "$CPPAD_DIR" ]; then
    echo "Could not find CppAD folder:"
    echo "$CPPAD_DIR"
    exit 1
fi
# -----------------------------------------------------------------------------
echo_eval rm -rf ${TMB_INCLUDE}/cppad
echo_eval cd ${CPPAD_DIR}/cppad; git clean -xdf
echo_eval rm -rf ${BUILD_CPPAD}
echo_eval mkdir ${BUILD_CPPAD}
echo_eval cd ${BUILD_CPPAD}; cmake -D eigen_prefix=${TMB_INCLUDE} -D cppad_testvector=eigen ..
echo_eval cp -r ${CPPAD_DIR}/cppad ${TMB_INCLUDE}
