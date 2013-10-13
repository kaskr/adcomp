rfiles = $(basename $(wildcard *.R))
cppfiles = $(basename $(wildcard *.cpp))
intersection= $(filter $(rfiles), $(cppfiles)) 
outputfiles = $(intersection:=.output.RData)

%.output.RData : %.R
	unset MAKEFLAGS; example=$(basename $<) R --vanilla < unittest.R

report: $(outputfiles)
	R --vanilla < unittest.R

clean :
	rm -f *.o *.so *~ core *.output.RData
