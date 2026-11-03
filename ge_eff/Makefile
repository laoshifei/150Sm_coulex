CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra $(shell root-config --cflags)
LDLIBS = $(shell root-config --libs) -lcrypto

all: source_spectra recover_legacy_spectra

source_spectra: source_spectra.cpp ../include/sm150_event_calibration.h ../include/sm150_parameter_paths.h
	$(CXX) $(CXXFLAGS) -I../include source_spectra.cpp -o $@ $(LDLIBS)

recover_legacy_spectra: recover_legacy_spectra.cpp ../include/sm150_event_calibration.h ../include/sm150_parameter_paths.h
	$(CXX) $(CXXFLAGS) -I../include recover_legacy_spectra.cpp -o $@ $(LDLIBS)

.PHONY: all clean
clean:
	rm -f source_spectra recover_legacy_spectra
