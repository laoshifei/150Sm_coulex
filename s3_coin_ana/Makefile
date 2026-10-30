out = ana

sourcefile = main.cpp analysis.h analysis.cpp s3_tree.h s3_tree.cpp set.h ../include/sm150_event_calibration.h

ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs) -lcrypto

CFLAGS = -Wall -Wextra -std=c++17 -O3 -ffp-contract=off -I../include -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

GXX = g++

all:$(out)

$(out): $(sourcefile) Makefile
	$(GXX) $(CFLAGS) $(ROOTCFLAGS) -o $@ $(filter %.cpp ,$(sourcefile)) $(ROOTLIBS)

clean:
	rm -f *~ *.o $(out)

.PHONY: all clean
