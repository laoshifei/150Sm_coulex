CXX = g++
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
CXXFLAGS = -Wall -Wextra -std=c++17 -O3 -ffp-contract=off -I../include

all: calibrate joint_calibrate

joint_calibrate: joint_calibrate.cpp ../include/sm150_event_calibration.h Makefile
	$(CXX) $(ROOTCFLAGS) $(CXXFLAGS) -o $@ joint_calibrate.cpp $(ROOTLIBS) -lcrypto

calibrate: calibrate.cpp ../include/sm150_event_calibration.h Makefile
	$(CXX) $(ROOTCFLAGS) $(CXXFLAGS) -o $@ calibrate.cpp $(ROOTLIBS) -lcrypto

.PHONY: all clean
clean:
	rm -f calibrate joint_calibrate
