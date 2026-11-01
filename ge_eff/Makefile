CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra $(shell root-config --cflags)
LDLIBS = $(shell root-config --libs) -lcrypto

source_spectra: source_spectra.cpp ../include/sm150_event_calibration.h
	$(CXX) $(CXXFLAGS) -I../include source_spectra.cpp -o $@ $(LDLIBS)

.PHONY: clean
clean:
	rm -f source_spectra
