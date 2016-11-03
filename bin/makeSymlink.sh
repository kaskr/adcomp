#/bin/bash

THISDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
TMBINST=`echo "cat(system.file(package='TMB'))" | R --slave`

if [ -d "${TMBINST}" ];
then
    rm -rf ${TMBINST}/include
    ln -s ${THISDIR}/../TMB/inst/include ${TMBINST}/include
    ls -l ${TMBINST}/include
else
    echo "TMB is not installed!"
fi
