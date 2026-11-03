CXX = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra $(shell root-config --cflags)
LDLIBS = $(shell root-config --libs) -lSpectrum -lMinuit2 -lcrypto

calibrate: main.cpp auto_cali.cpp get_par.cpp compute_average.cpp joint_source.h ../include/sm150_source_calibration.h ../include/sm150_event_calibration.h ../include/sm150_parameter_paths.h
	$(CXX) $(CXXFLAGS) main.cpp -o $@ $(LDLIBS)

.PHONY: clean
clean:
	rm -f calibrate
