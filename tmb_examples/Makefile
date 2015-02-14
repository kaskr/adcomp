excludes = lr_test spde_aniso spde
rfiles = $(basename $(wildcard *.R))
cppfiles = $(basename $(wildcard *.cpp))
all_examples = $(filter $(rfiles), $(cppfiles))
examples = $(filter-out $(excludes), $(all_examples))
outputfiles = $(examples:=.output.RData)
profiletargets = $(examples:=.profile)

%.output.RData : %.R %.cpp
	unset MAKEFLAGS; example=$(basename $<) R --vanilla < unittest.R

all: $(outputfiles)
	make report

$(all_examples) : % : %.output.RData
	

report:
	R --vanilla < unittest.R

clean :
	rm -rf *.o *.so *~ core *.output.RData Rplots.pdf *.profile

%.profile : %.R %.cpp
	example=$(basename $<) R --vanilla < profiler.R

profile_all: $(outputfiles) $(profiletargets)

