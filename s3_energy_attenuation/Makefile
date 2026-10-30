CXX = g++
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
CXXFLAGS = -Wall -Wextra -std=c++17 -O3 -ffp-contract=off -I../include

all: extract calibrate conditioned_extract conditioned_calibrate

extract: tree2hist.cpp ../include/sm150_event_calibration.h Makefile
	$(CXX) $(ROOTCFLAGS) $(CXXFLAGS) -o $@ tree2hist.cpp $(ROOTLIBS) -lcrypto

.PHONY: all clean
calibrate: fit_xy_new2.cpp ../include/sm150_event_calibration.h Makefile
	$(CXX) $(ROOTCFLAGS) $(CXXFLAGS) -o $@ fit_xy_new2.cpp $(ROOTLIBS) -lcrypto

clean:
	rm -f extract calibrate conditioned_extract conditioned_calibrate

conditioned_extract: conditioned_extract.cpp ../include/sm150_event_calibration.h Makefile
	$(CXX) $(ROOTCFLAGS) $(CXXFLAGS) -o $@ conditioned_extract.cpp $(ROOTLIBS) -lcrypto

conditioned_calibrate: conditioned_calibrate.cpp ../include/sm150_event_calibration.h Makefile
	$(CXX) $(ROOTCFLAGS) $(CXXFLAGS) -o $@ conditioned_calibrate.cpp $(ROOTLIBS) -lcrypto
