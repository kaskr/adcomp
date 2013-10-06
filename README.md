adcomp
======
* Development of RcppAD package.
* AD computation with RcppAD.

To install package
==================
* make install

To build manual
===============
* make pdf

Alternative platforms
=====================

Windows
-------
Tested to work on 64 bit R with latest [Rtools](http://cran.r-project.org/bin/windows/Rtools/).
Currently only the non-parallel examples run.

Mac OS X
--------
Tested to work with both llvm-gcc-4.2 and clang. Fortran compiler [libraries](http://cran.r-project.org/bin/macosx/tools) must be installed. According to [R admin manual](http://www.cran.r-project.org/doc/manuals/R-admin.html#OS-X) "the OpenMP support in this version of gcc is problematic, and the alternative, clang, has no OpenMP support". So, parallel templates with OS X will require some work.
