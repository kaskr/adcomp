examples=nmix mvrw spatial sdv_multi thetalog orange_big
outputfiles = $(examples:=.output.RData)
all: $(outputfiles)
%.output.RData :
	example=$(basename $@) R --vanilla < run.R

