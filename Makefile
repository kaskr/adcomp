PACKAGE=TMB
VERSION=1.1
TARBALL=${PACKAGE}_${VERSION}.tar.gz
ZIPFILE=${PACKAGE}_${VERSION}.zip

all:
	make doc-update
	make build-package
	make install
	make pdf

doc-update:
	echo "library(roxygen2);roxygenize(\"$(PACKAGE)\",roclets = c(\"collate\", \"rd\"))" | R --slave

build-package:
	R CMD build --resave-data=no $(PACKAGE)

install:
	make build-package
	R CMD INSTALL --preclean $(TARBALL)

install-metis:
	make build-package
	LIBCHOLMOD=/usr/lib/libcholmod.so.3.1.0 R CMD INSTALL --preclean $(TARBALL)


unexport TEXINPUTS
pdf:
	rm -f $(PACKAGE).pdf
	R CMD Rd2pdf --no-preview $(PACKAGE)

check:
	R CMD check $(PACKAGE)

unlock:
	rm -rf ${R_LIBS}/00LOCK-TMB

.PHONY: dox
dox:
	cd dox; doxygen

dox-clean:
	cd dox; rm -rf html latex

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

install-metis-full: $(SUITESPARSE).tar.gz $(METIS).tar.gz
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
	make build-package
	LIBCHOLMOD=`pwd`/SuiteSparse/libcholmod.so R CMD INSTALL --preclean $(TARBALL)
