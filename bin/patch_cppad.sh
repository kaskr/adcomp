#! /bin/bash -e
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
    echo $*
    eval $*
}
# -----------------------------------------------------------------------------
THISDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
TMBDIR=${THISDIR}/..
OPTIONS="--batch"
PATCH="echo_eval patch ${OPTIONS} -d${TMBDIR} -p1 < "
# -----------------------------------------------------------------------------
# Go to folder with patch dir
echo_eval cd ${THISDIR}
# -----------------------------------------------------------------------------
# Apply patches
${PATCH} patch/CTOR-workarounds.patch
${PATCH} patch/CppAD-jacobian-prefer-reverse-mode-over-forward-mode.patch
${PATCH} patch/Parallel-checkpoint-still-needed.patch
