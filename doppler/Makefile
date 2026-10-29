out = doppler

sourcefile = main.cpp doppler.h doppler.cpp event.h event.cpp angle_info.h angle_info.cpp set.h

ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs) -lSpectrum -lcrypto

CFLAGS = -g -Wall -Wextra -std=c++17 -O3 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

GXX = g++

all:$(out)

$(out): $(sourcefile)
	$(GXX) $(CFLAGS) $(ROOTCFLAGS) -o $@ $(filter %.cpp ,$(sourcefile)) $(ROOTLIBS)

clean:
	rm -f *~ *.o $(out)
