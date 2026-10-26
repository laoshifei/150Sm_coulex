CXX = g++
CXXFLAGS = -O2 -g -Wall -Wextra -std=c++17 -I../include
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs) -lcrypto
.PHONY: all clean
all: build fit
build: main.cpp build.cpp build.h raw.cpp raw.h ../include/sm150_event_calibration.h ../include/sm150_time_alignment.h Makefile
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) main.cpp build.cpp raw.cpp -o $@.tmp $(ROOTLIBS)
	mv -f $@.tmp $@
fit: fit_time_walk.cpp candidate.h ../include/sm150_event_calibration.h ../include/sm150_time_alignment.h Makefile
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) fit_time_walk.cpp -o $@.tmp $(ROOTLIBS) -lMinuit2
	mv -f $@.tmp $@
clean:
	rm -f build fit *.o
