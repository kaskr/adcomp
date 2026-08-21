tarballs = $(wildcard *.tar.gz)
pkgs = $(basename $(basename $(tarballs)))
Rchecks = $(pkgs:=.Rcheck)

all: $(Rchecks)

## Get tarballs:
init:
	R --slave < init.R

%.Rcheck : %.tar.gz
	-R CMD check $<

clean:
	rm -rf *.Rcheck

## Inspect:
status:
	grep Status */00check.log
can-be-installed:
	grep 'can be installed' *.Rcheck/00check*
