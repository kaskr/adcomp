rfiles = $(basename $(wildcard *.R))
cppfiles = $(basename $(wildcard *.cpp))
intersection= $(filter $(rfiles), $(cppfiles)) 
outputfiles = $(intersection:=.output.RData)
profiletargets = $(intersection:=.profile)

%.output.RData : %.R %.cpp
	unset MAKEFLAGS; example=$(basename $<) R --vanilla < unittest.R

all: $(outputfiles)
	make report

$(intersection) : % : %.output.RData
	

report:
	R --vanilla < unittest.R

clean :
	rm -rf *.o *.so *~ core *.output.RData Rplots.pdf *.profile

%.profile : %.R %.cpp
	example=$(basename $<) R --vanilla < profiler.R

profile_all: $(outputfiles) $(profiletargets)

