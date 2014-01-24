examples=nmix mvrw spatial sdv_multi
outputfiles = $(examples:=.output.RData)
all: $(outputfiles)
%.output.RData :
	example=$(basename $@) R --vanilla < run.R

