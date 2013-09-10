.SUFFIXES:
.SUFFIXES: .cpp

rfiles = $(basename $(wildcard *.R))
cppfiles = $(basename $(wildcard *.cpp))

.cpp:
	echo "library(RcppAD);compile(\"$@.cpp\")" | R --vanilla

$(foreach doc,$(cppfiles),$(eval $(doc):$(doc).cpp))

all: $(cppfiles)

run: all
	$(foreach doc,$(rfiles),R --vanilla < $(doc).R > $(doc).out;)
	rm -f examples.out
	cat *.out > examples.out

clean :
	rm -f *.o *.so *.out *~ core *.RData




bug:
	g++ -I/appl/R/R-2.15.2/lib64/R/include -DNDEBUG  -I../include    -fpic  -g -O0  -c randomregression.cpp -o randomregression.o
	g++ -shared -L/usr/local/lib64 -o randomregression.so randomregression.o -L/appl/R/R-2.15.2/lib64/R/lib -lR
