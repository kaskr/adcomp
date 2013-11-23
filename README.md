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

Linux specific notes
====================

Metis orderings
---------------
For large 3D random field models the ordering algorithms shipping with R's Matrix package are far from optimal. To get better orderings available run the following in the terminal:

* sudo apt-get install libsuitesparse-metis-3.1.0

This will install a more complete version of CHOLMOD with more orderings available. Then install TMB like this:

* make install-metis

For a quick example of how to use it start R, load TMB and run

* runExample("ar1xar1")

Issues with library unloading
-----------------------------
On recent versions of gcc the following problem may be encountered: When the user cpp file is changed, re-compiled and re-loaded, the changes does not take place. To see if you are affected by this issue, assuming your compiled DLL is called "mymodel.so", try running:

* readelf -s mymodel.so | grep UNIQUE

If this gives a lot of output it is not possible to unload the library, and R will have to be restarted every time the model is re-compiled.
A workaround is to use clang++ instead of gcc.

