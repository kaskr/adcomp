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
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\",roclets = c(\"collate\", \"rd\"))" | $(R) --slave

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
## Select version that match R's Matrix package
SUITESPARSE = SuiteSparse-4.2.1
METIS = metis-4.0.3
WGET = curl -O
OS = $(shell uname)

$(SUITESPARSE).tar.gz :
	$(WGET) http://faculty.cse.tamu.edu/davis/SuiteSparse/$(SUITESPARSE).tar.gz

$(METIS).tar.gz :
	$(WGET) http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/$(METIS).tar.gz

SuiteSparse: $(SUITESPARSE).tar.gz $(METIS).tar.gz
	tar zxfv $(SUITESPARSE).tar.gz
	cd SuiteSparse; cp ../$(METIS).tar.gz .
	cd SuiteSparse; tar zxfv $(METIS).tar.gz
	cd SuiteSparse; ln -s $(METIS) metis-4.0
## Edit "metis-4.0/Makefile.in" with COPTIONS = -fPIC
	cd SuiteSparse; sed -i.backup s/COPTIONS\ =/COPTIONS\ =\ -fPIC/g metis-4.0/Makefile.in
	cd SuiteSparse; cd metis-4.0 && make
	cd SuiteSparse; make library
## Restore object files so we can make .so
	cd SuiteSparse; cd SuiteSparse_config; ar vx *.a
	cd SuiteSparse; cd CCOLAMD/Lib; ar vx *.a
	cd SuiteSparse; cd COLAMD/Lib; ar vx *.a
	cd SuiteSparse; gcc -shared -o libcholmod.so SuiteSparse_config/SuiteSparse_config.o CHOLMOD/Lib/*.o AMD/Lib/*.o CAMD/Lib/*.o CCOLAMD/Lib/*.o COLAMD/Lib/*.o metis-4.0/Lib/*.o `R CMD config BLAS_LIBS` `R CMD config LAPACK_LIBS`
	if [ "$(OS)" = "Darwin" ]; then						\
		cd SuiteSparse;							\
		install_name_tool -id `pwd`/libcholmod.so libcholmod.so;	\
	fi

install-metis-full: SuiteSparse
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
	echo ".onAttach <- function(lib, pkg) {"                                                      >> TMB/R/zzz.R
	echo "  exfolder <- system.file(\"examples\", package = \"TMB\")"                             >> TMB/R/zzz.R
	echo "  dll <- paste0(exfolder, Sys.getenv(\"R_ARCH\"), \"/simple\", .Platform\$$dynlib.ext)" >> TMB/R/zzz.R
	echo "  if(!file.exists(dll)) runExample(\"simple\", dontrun=TRUE)"                           >> TMB/R/zzz.R
	echo "}"                                                                                      >> TMB/R/zzz.R
	make eliminate-cout
	make doc-update
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
