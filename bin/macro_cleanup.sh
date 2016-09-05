#!/bin/bash
# Input:  C++ file
# Output: File 'undefs.h' that cleans up macros specific for input file.
if [ "$1" == "" ]
then
	echo "Must specify a filename"
	exit 1
fi

echo "#include <R.h>"          >  tmp_prog1.cpp
echo "#include <Rinternals.h>" >> tmp_prog1.cpp

echo "#include <R.h>"          >  tmp_prog2.cpp
echo "#include <Rinternals.h>" >> tmp_prog2.cpp
echo "#include \"$1\""         >> tmp_prog2.cpp

gcc -I/usr/share/R/include -E -dM tmp_prog1.cpp > out1
gcc -I/usr/share/R/include -E -dM tmp_prog2.cpp > out2

diff out1 out2 > out.diff

grep '^<' out.diff | sed 's/^</!!! RE-DEFINED !!!: /g'
grep '^>' out.diff | sed 's/^.*define //g' | sed 's/[ (].*//g' | sed 's/^/#undef /g' | grep -v TINY_AD > undefs.h
echo "Generated 'undefs.h'"
