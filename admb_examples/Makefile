examples=nmix mvrw spatial thetalog socatt sdv_multi orange_big sam
outputfiles = $(examples:=.output.RData) longlinreg.output.RData
all: $(outputfiles)
%.output.RData :
	example=$(basename $@) R --vanilla < run.R

longlinreg.output.RData:
	R --vanilla < run_longlinreg.R
