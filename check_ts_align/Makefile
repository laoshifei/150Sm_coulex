out = build

sourcefile = main.cpp build.h build.cpp set.h ../include/sm150_time_alignment.h

ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs) -lSpectrum

CFLAGS = -g -Wall -std=c++11 -O2 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

GXX = g++

.PHONY: all legacy clean

all: check_build

# The original raw-data builder is retained only for historical reproduction.
legacy: $(out)

check_build: check_build.cpp Makefile
	$(GXX) $(ROOTCFLAGS) -Wall -Wextra -std=c++17 -O3 -o $@ check_build.cpp $(ROOTLIBS) -lcrypto

$(out): $(sourcefile)
	$(GXX) $(CFLAGS) $(ROOTCFLAGS) -o $@ $(filter %.cpp ,$(sourcefile)) $(ROOTLIBS)

clean:
	rm -f *~ *.o $(out) check_build
