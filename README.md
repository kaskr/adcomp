Template Model Builder (TMB)
============================
TMB is an R package with functionality similar to ADMB.
It requires R at least version 3.0.0 and development tools needed to install R packages from source.
The package is installed from the command line by entering the adcomb folder and typing

* make install

To build the manual type

* make pdf

Once the package is successfully installed it can be tested by starting R and running

    library(TMB)
    runExample(all=TRUE)

Alternative platforms
=====================

Windows
-------
Tested to work on 64 bit R with latest [Rtools](http://cran.r-project.org/bin/windows/Rtools/).
All examples run (including the parallel), however it may be needed to manually edit the PKG_LIBS of the Makevars file in order to correctly link to the OpenMP runtime library. Otherwise it will not be possible to change the number of threads from within R. Note that the folder where the package is installed is not allowed to contain spaces.

Mac OS X
--------
Tested to work with both llvm-gcc-4.2 and clang. Fortran compiler [libraries](http://cran.r-project.org/bin/macosx/tools) must be installed. According to [R admin manual](http://www.cran.r-project.org/doc/manuals/R-admin.html#OS-X) "the OpenMP support in this version of gcc is problematic, and the alternative, clang, has no OpenMP support". So, parallel templates with OS X will require a different compiler installed.
