examples=nmix mvrw spatial thetalog socatt sdv_multi orange_big
outputfiles = $(examples:=.output.RData)
all: $(outputfiles)
%.output.RData :
	example=$(basename $@) R --vanilla < run.R

