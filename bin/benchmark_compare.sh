#!/bin/bash
# -----------------------------------------------------------------------------
# bash function that rename files by regexp. Ignoring 1st arg.
rename() {
    for f in ${@:3}; do mv "$f" $(echo "$f" | sed $2); done
}
# ----------------------------------------------------------------------------

BRANCH1=$1
BRANCH2=$2
THISDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BENCHMARK=${THISDIR}/../benchmark
BENCHMARK_DIR1=${BENCHMARK}/${BRANCH1}
BENCHMARK_DIR2=${BENCHMARK}/${BRANCH2}
BENCHMARK_COMPARE=${BENCHMARK}/compare-${BRANCH1}-${BRANCH2}
VALID_ARGS=`ls ${BENCHMARK} | grep -v "^compare-"`

## Check inputs
if [[ ! $VALID_ARGS =~ (^|[[:space:]])"$1"($|[[:space:]]) ]] ; then
    echo "1st argument must be one of:" $VALID_ARGS;
    exit 1;
fi
if [[ ! $VALID_ARGS =~ (^|[[:space:]])"$2"($|[[:space:]]) ]] ; then
    echo "2nd argument must be one of:" $VALID_ARGS;
    exit 1;
fi
## OK

mkdir -p ${BENCHMARK_COMPARE}
rm -rf ${BENCHMARK_COMPARE}/*

## Copy branch1 and branch2 results to common folder
## Setup tools
cp ${BENCHMARK_DIR1}/Makefile ${BENCHMARK_COMPARE}
cd ${BENCHMARK_DIR1}
find . -name '*.cpp'          | cpio -pdm ${BENCHMARK_COMPARE}
find . -name '*.R'            | cpio -pdm ${BENCHMARK_COMPARE}
find . -name '*.output.RData' | cpio -pdm ${BENCHMARK_COMPARE}
find . -name '*.logpid'       | cpio -pdm ${BENCHMARK_COMPARE}

## - branch1 becomes 'expected' results
cd ${BENCHMARK_COMPARE}
rename -f 's/\.output.RData$/\.expected.RData/' *.output.RData */*.output.RData
rename -f 's/\.logpid$/\.logpid1/' *.logpid */*.logpid

## - branch2 becomes 'output' results
cd ${BENCHMARK_DIR2}
find . -name '*.output.RData' | cpio -pdm ${BENCHMARK_COMPARE}
find . -name '*.logpid'       | cpio -pdm ${BENCHMARK_COMPARE}
cd ${BENCHMARK_COMPARE}
rename -f 's/\.logpid$/\.logpid2/' *.logpid */*.logpid

## Table comparision (REPORT.md)
cd ${BENCHMARK_COMPARE}
make report

## Fig memory profile (memory.pdf)
echo '
files1 <- dir(pattern="*.logpid1$", recursive=TRUE)
files2 <- dir(pattern="*.logpid2$", recursive=TRUE)
files  <- intersect(sub(".logpid1$","",files1),
                    sub(".logpid2$","",files2))
files1 <- paste0(files, ".logpid1")
files2 <- paste0(files, ".logpid2")
f <- function(f1, f2) {
    x <- read.table(f1, TRUE)
    y <- read.table(f2, TRUE)
    m <- matrix(NA, max(nrow(x), nrow(y)), 2)
    m[seq_len(nrow(x)),1] <- x$pmem
    m[seq_len(nrow(y)),2] <- y$pmem
    t <- seq_len(nrow(m)) - 1
    t <- t / 10
    matplot(t, m, type="b", ylab="Memory (percent)", xlab="Time (sec)")
    legend("bottomright", legend=c("BRANCH1", "BRANCH2"), col=1:2, lty=1, lwd=2)
    title(sub(".logpid1$", "", f1))
}
pdf("memory.pdf")
invisible(Map(f, files1, files2))
dev.off()
' |
sed s/BRANCH1/${BRANCH1}/g |
sed s/BRANCH2/${BRANCH2}/g |
R --slave

echo
echo "Output generated:"
echo ${BENCHMARK_COMPARE}/REPORT.md
echo ${BENCHMARK_COMPARE}/memory.pdf
