PACKAGE=TMB
VERSION := $(shell sed -n '/^Version: /s///p' TMB/DESCRIPTION)
DATE := $(shell sed -n '/^Date: /s///p' TMB/DESCRIPTION)
TARBALL=${PACKAGE}_${VERSION}.tar.gz
ZIPFILE=${PACKAGE}_${VERSION}.zip

# Allow e.g. "make R=R-devel install"
R=R

all:
	make doc-update
	make build-package
	make install
	make pdf

doc-update:
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\",roclets = c(\"collate\", \"rd\"), load_code=load_source)" | $(R) --slave

build-package:
	$(R) CMD build --resave-data=no $(PACKAGE)

install:
	make build-package
	$(R) CMD INSTALL --preclean $(TARBALL)

install-metis:
	make build-package
	LIBCHOLMOD=/usr/lib/libcholmod.so.3.1.0 $(R) CMD INSTALL --preclean $(TARBALL)


unexport TEXINPUTS
pdf:
	rm -f $(PACKAGE).pdf
	$(R) CMD Rd2pdf --no-preview $(PACKAGE)

check:
	$(R) CMD check $(PACKAGE)

unlock:
	rm -rf `Rscript --vanilla -e 'writeLines(.Library)'`/00LOCK-TMB

## Alternative 'install-metis': Get source code and build...
## Select version matching R's Matrix::.SuiteSparse_version()
## SUITESPARSE_VERSION = 4.2.1
SUITESPARSE_VERSION = 5.7.1
## METIS = metis-4.0.3
WGET = wget
OS = $(shell uname)

## Download SuiteSparse
v$(SUITESPARSE_VERSION).tar.gz :
	$(WGET) https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/v$(SUITESPARSE_VERSION).tar.gz

## Unpack SuiteSparse
SuiteSparse: v$(SUITESPARSE_VERSION).tar.gz
	tar zxfv v$(SUITESPARSE_VERSION).tar.gz
	ln -sf SuiteSparse-$(SUITESPARSE_VERSION) SuiteSparse

##$(METIS).tar.gz :
##	$(WGET) http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/$(METIS).tar.gz

## Build METIS and CHOLMOD
SuiteSparse/libcholmod.so: SuiteSparse
	( cd SuiteSparse && $(MAKE) metis)
	( cd SuiteSparse/SuiteSparse_config && $(MAKE) )
	( cd SuiteSparse/CAMD && $(MAKE) )
	( cd SuiteSparse/AMD && $(MAKE) )
	( cd SuiteSparse/CCOLAMD && $(MAKE) )
	( cd SuiteSparse/COLAMD && $(MAKE) )
	( cd SuiteSparse/CHOLMOD && BLAS=`R CMD config BLAS_LIBS` $(MAKE) library )
	( cd SuiteSparse &&  gcc -shared -o libcholmod.so `find . -type f -name "*.o" | grep -v programs` )
	if [ "$(OS)" = "Darwin" ]; then						\
		cd SuiteSparse;							\
		install_name_tool -id `pwd`/libcholmod.so libcholmod.so;	\
	fi

install-metis-full: SuiteSparse/libcholmod.so
	make build-package
	LIBCHOLMOD=`pwd`/SuiteSparse/libcholmod.so R CMD INSTALL --preclean $(TARBALL)

## Get a rough changelog since most recent github revision tag
## (Use as starting point when updating NEWS file)
## NOTE: Run *after* updating version and date in DESCRIPTION.
changelog:
	echo; \
	echo "------------------------------------------------------------------------"; \
	echo TMB $(VERSION) \($(DATE)\); \
	echo "------------------------------------------------------------------------"; \
	echo; \
	git --no-pager log --format="o %B" `git describe --abbrev=0 --tags`..HEAD | sed s/^-/\ \ -/g

## The CRAN version must be customized a bit...
## FIXME: Is it possible to get 'Makevars' POSIX compliant without
## losing the current flexibility e.g. 'make install-metis'? (Main
## obstacle is lack of 'ifdef' in POSIX make).
## FIXME: 'LinkingTo RcppEigen' is not really right but perhaps the
## best we can do to assert RcppEigen is installed? (we need the
## Eigen headers, nothing else)
eliminate-cout:
	cd TMB/inst/include; sed -i /.*using\ std::cout.*/d cppad/*.hpp cppad/*/*.hpp
	cd TMB/inst/include; sed -i s/[std:]*cout/Rcout/g cppad/*.hpp cppad/*/*.hpp
	cd TMB/inst/include; sed -i s/std..cout/Rcout/g ./*.hpp tmbutils/*.hpp
	git checkout TMB/inst/include/Rstream.hpp
cran-version:
	cd TMB; git clean -xdf
	sed -i 's/^LinkingTo.*/LinkingTo: Matrix, RcppEigen/' TMB/DESCRIPTION
	rm -rf TMB/inst/include/Eigen
	rm -rf TMB/inst/include/unsupported
	echo "PKG_LIBS = \$$(LAPACK_LIBS) \$$(BLAS_LIBS) \$$(FLIBS) \$$(SHLIB_OPENMP_CFLAGS)" > TMB/src/Makevars
	echo "PKG_CFLAGS = \$$(SHLIB_OPENMP_CFLAGS)"                                         >> TMB/src/Makevars
	sed -i /^SystemRequirements.*/d TMB/DESCRIPTION
	echo ".onAttach <- function(lib, pkg) {"                                                          >> TMB/R/zzz.R
	echo "  exfolder <- system.file(\"examples\", package = \"TMB\")"                                 >> TMB/R/zzz.R
	echo "  dll <- paste0(exfolder, Sys.getenv(\"R_ARCH\"), \"/simple\", .Platform\$$dynlib.ext)"     >> TMB/R/zzz.R
	echo "  if(!file.exists(dll)) runExample(\"simple\", dontrun=TRUE, eigen.disable.warnings=FALSE)" >> TMB/R/zzz.R
	echo "}"                                                                                          >> TMB/R/zzz.R
	make eliminate-cout
	make build-package

##########################################################
## For travis tests
test-tmb_syntax:
	$(R) --version
	cd tmb_syntax; make test

test-tmb_examples:
	$(R) --version
	cd tmb_examples; make test

doxygen:
	cd dox; make all

cran-check:
	$(R) CMD check --as-cran $(TARBALL)
