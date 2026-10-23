out = build

sourcefile = main.cpp build.h build.cpp ../include/sm150_time_alignment.h Makefile

ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs) -lSpectrum

CFLAGS = -Wall -std=c++17 -O3 -ffp-contract=off -I../include -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

GXX = g++

all:$(out)

$(out): $(sourcefile)
	$(GXX) $(ROOTCFLAGS) $(CFLAGS) -o $@ $(filter %.cpp ,$(sourcefile)) $(ROOTLIBS) -lcrypto

clean:
	rm -f *~ *.o $(out)
