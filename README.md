Template Model Builder (TMB)
============================
TMB is an R package with functionality similar to ADMB.
It requires R at least version 3.0.0 and development tools needed to install R packages from source.
It is recommended to install TMB into your local R package library (if you do not yet have a local R package library, create one by running ```install.packages("")``` from the R prompt and follow the instructions).
The package is installed from the command line by entering the adcomp folder and typing

* make install

To build the user manual type

* make pdf

Once the package is successfully installed it can be tested by starting R and running

    library(TMB)
    runExample(all=TRUE)

To build API-Function reference (requires that you have "doxygen" installed) type

* make dox

Once the documentation is buildt open dox/html/index.html in your web browser. 
Use the search field to find functions and their documentation.

It is recommended to test that models can be changed, re-compiled and re-loaded without problems:

    setwd("tmb_syntax")
    source("test_reload.R")

If everything works as intended this script should display "TRUE".
The script has only been observed to fail for certain combinations of Linux and the gcc compiler, see below.

Compiler warnings
-----------------
Some R configurations have the '-Wall' compiler flag enabled by default resulting in a large number of useless warnings. You can list your C++ compiler flags from the terminal with ```R CMD config CXXFLAGS```.
It is recommended to remove '-Wall' from this list, e.g. by running the following in a terminal:

    mkdir -p ~/.R
    echo CXXFLAGS=`R CMD config CXXFLAGS | sed s/-Wall//g` >> ~/.R/Makevars


Alternative platforms
=====================

Windows
-------
Tested to work on 64 bit R with latest [Rtools](http://cran.r-project.org/bin/windows/Rtools/). Currently not working with 32 bit R.

_Install instructions_

1. Start 64 bit R and change working directory to the (cloned or unzipped) ```adcomp``` folder.

2. From R run: ```source("install_windows.R")```

The required Rtools will be downloaded and installed. Note that the PATH variable need not be changed by the installer or otherwise edited. The PATH will be automatically set for each TMB session.

_Status_

- Parallel user templates work, including changing the number of threads from R.
- Filenames and folders with spaces should be ok.
- -Wall flag disabled by default.

Mac OS X
--------
Tested to work with both llvm-gcc-4.2 and clang. Fortran compiler [libraries](http://cran.r-project.org/bin/macosx/tools) must be installed. According to [R admin manual](http://www.cran.r-project.org/doc/manuals/R-admin.html#OS-X) "the OpenMP support in this version of gcc is problematic, and the alternative, clang, has no OpenMP support". So, parallel templates with OS X will require a different compiler installed.

On Mavericks (Mac OS 10.9) the Fortran compiler on the CRAN website does not work. More information and the solution can be found [here](http://www.thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error). The easiest solution is to install the appropriate fortran libraries from r.research.att via the command line in Terminal:

```
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

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

