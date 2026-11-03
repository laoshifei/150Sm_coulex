out = sort

sourcefile = main.cpp sort.h sort.cpp raw.h raw.cpp set.h ../include/sm150_time_alignment.h ../include/sm150_parameter_paths.h

ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs) -lSpectrum -lcrypto

CFLAGS = -g -Wall -std=c++17 -O2 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

GXX = g++

all:$(out)

$(out): $(sourcefile)
	$(GXX) $(CFLAGS) $(ROOTCFLAGS) -o $@ $(filter %.cpp ,$(sourcefile)) $(ROOTLIBS)

clean:
	rm -f *~ *.o $(out)
