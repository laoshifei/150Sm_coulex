CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -std=c++17 -ffp-contract=off $(shell root-config --cflags)
LDLIBS = $(shell root-config --libs) -lcrypto

histogram: main.cpp ../include/sm150_event_calibration.h Makefile
	$(CXX) $(CXXFLAGS) -I../include main.cpp -o $@ $(LDLIBS)

.PHONY: clean
clean:
	rm -f histogram
