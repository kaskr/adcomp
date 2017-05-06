excludes = lr_test spde_aniso spde
rfiles = $(basename $(wildcard *.R)) $(basename $(wildcard */*.R))
cppfiles = $(basename $(wildcard *.cpp)) $(basename $(wildcard */*.cpp))
all_examples = $(filter $(rfiles), $(cppfiles))
examples = $(filter-out $(excludes), $(all_examples))
outputfiles = $(examples:=.output.RData)
profiletargets = $(examples:=.profile)

%.output.RData : %.R %.cpp
	unset MAKEFLAGS; example=$(basename $<) R --slave < tools/unittest.R

all: $(outputfiles)
	make report

$(all_examples) : % : %.output.RData
	

report:
	report_level=0 R --slave < tools/unittest.R

report-full:
	report_level=1 R --slave < tools/unittest.R

clean :
	rm -rf *.o *.so *.dll *~ core *.output.RData *.pdf *.profile
	rm -rf */*.o */*.so */*.dll */*~ core */*.output.RData */*.pdf */*.profile

%.profile : %.R %.cpp
	example=$(basename $<) R --slave < tools/profiler.R

%.logpid : %.R %.cpp
	example=$(basename $<) R --slave < tools/logpid.R

profile_all: $(outputfiles) $(profiletargets)

##########################################################
## For travis tests
test:
	make $(examples) | grep slave
	make report
	cat REPORT.md
	if grep -q FALSE REPORT.md; then exit 1; fi;
	if grep -q NA REPORT.md; then exit 1; fi;
