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
Not yet fully tested. Main difference from linux builds is that the mac toolset lacks Fortran libraries (used to link to BLAS) and use the older gcc-4.2.
