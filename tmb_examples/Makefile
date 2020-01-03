excludes = lr_test
rfiles = $(basename $(wildcard *.R)) $(basename $(wildcard */*.R))
cppfiles = $(basename $(wildcard *.cpp)) $(basename $(wildcard */*.cpp))
all_examples = $(filter $(rfiles), $(cppfiles))
examples = $(filter-out $(excludes), $(all_examples))
outputfiles = $(examples:=.output.RData)
profiletargets = $(examples:=.profile)
logpidtargets = $(examples:=.logpid)

# Allow e.g. "make R=R-devel"
R=R

%.output.RData : %.R %.cpp
	unset MAKEFLAGS; example=$(basename $<) $(R) --slave < tools/unittest.R

all: $(outputfiles)
	make report

$(all_examples) : % : %.output.RData
	

report:
	report_level=0 $(R) --slave < tools/unittest.R

report-full:
	report_level=1 $(R) --slave < tools/unittest.R

clean :
	rm -rf *.o *.so *.dll *~ core *.output.RData *.pdf *.profile *.memprofile
	rm -rf */*.o */*.so */*.dll */*~ core */*.output.RData */*.pdf */*.profile */*.memprofile

%.profile : %.R %.cpp
	example=$(basename $<) type=hotspots $(R) --slave < tools/profiler.R

%.memprofile : %.R %.cpp
	example=$(basename $<) type=mi1 $(R) --slave < tools/profiler.R

%.logpid : %.R %.cpp
	example=$(basename $<) $(R) --slave < tools/logpid.R

profile_all: $(outputfiles) $(profiletargets)

logpid_all: $(outputfiles) $(logpidtargets)

##########################################################
## For travis tests
test:
	make $(examples) | grep slave
	make report
	cat REPORT.md
	if grep -q FALSE REPORT.md; then exit 1; fi;
	if grep -q NA REPORT.md; then exit 1; fi;
