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

