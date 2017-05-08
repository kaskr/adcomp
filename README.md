Template Model Builder (TMB)
============================

[![Build Status](https://travis-ci.org/kaskr/adcomp.svg?branch=master)](https://travis-ci.org/kaskr/adcomp)

TMB is an R package with functionality similar to ADMB.
It requires R at least version 3.0.0 and development tools needed to install R packages from source.
Standard install instructions are available [here](https://github.com/kaskr/adcomp/wiki/Download).

# Install the development version (Linux)

It is recommended to install TMB into your local R package library (if you do not yet have a local R package library, create one by running ```install.packages("")``` from the R prompt and follow the instructions).
The package is installed from the command line by entering the adcomp folder and typing

```shell
make install
```

To build the user manual type

```shell
make pdf
```

Once the package is successfully installed it can be tested by starting R and running

```R
library(TMB)
runExample(all=TRUE)
```

To build API-Function reference see instructions [here](dox/README.md).
Once the documentation is built open `dox/html/index.html` in your web browser.
Use the search field to find functions and their documentation.

It is recommended to test that models can be changed, re-compiled and re-loaded without problems:

```R
setwd("tmb_syntax")
source("test_reload.R")
```

If everything works as intended this script should display "TRUE".
The script has only been observed to fail for certain combinations of Linux and the gcc compiler, see below.

# Advanced notes

## Metis orderings

For large 3D random field models the ordering algorithms shipping with R's Matrix package are far from optimal. To get better orderings available run the extended TMB install script (Linux/Mac only):

```shell
make install-metis-full
```

The new orderings are available from R through the function `runSymbolicAnalysis()`.
For a quick example of how to use it try the [ar1_4D](./tmb_examples/ar1_4D.R) example:

```shell
cd tmb_examples
make ar1_4D
```

You should see a significant speedup.

## Issues with library unloading

On recent versions of gcc the following problem may be encountered: When the user cpp file is changed, re-compiled and re-loaded, the changes do not take place. To see if you are affected by this issue, assuming your compiled DLL is called "mymodel.so", try running:

```shell
readelf -s mymodel.so | grep UNIQUE
```

If this gives any output it is not possible to unload the library, and R will have to be restarted every time the model is re-compiled.
There are at least two alternative solutions to this problem:

1. Use gcc with compilation flag ```-fno-gnu-unique``` (version 4.8.3 and newer): Add ```CXX = g++ -fno-gnu-unique``` to a file ```~/.R/Makevars``` (create it if it doesn't exist).

2. Use the clang compiler instead of gcc: Install clang and add ```CXX = clang++``` to a file ```~/.R/Makevars``` (create it if it doesn't exist).

Note: If you have precompiled TMB using ```precompile()``` prior to (1) or (2) you must repeat the precompilation with the new compiler settings.
