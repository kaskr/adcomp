examples=nmix mvrw
outputfiles = $(examples:=.output.RData)
all: $(outputfiles)
%.output.RData :
	example=$(basename $@) R --vanilla < run.R

