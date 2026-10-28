CXX = g++
CXXFLAGS = -g -Wall -Wextra -O2 -std=c++17 -I../include
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs) -lcrypto
.PHONY: all clean
all: analysis
analysis: main.cpp analysis.cpp analysis.h raw.cpp raw.h set.h ../include/sm150_time_alignment.h ../include/sm150_event_calibration.h
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) main.cpp analysis.cpp raw.cpp -o $@ $(ROOTLIBS)
clean:
	rm -f analysis *.o
