out = doppler

sourcefile = main.cpp doppler.h doppler.cpp event.h event.cpp angle_info.h angle_info.cpp set.h ../include/sm150_parameter_paths.h

ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs) -lSpectrum -lcrypto

CFLAGS = -g -Wall -Wextra -std=c++17 -O3 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

GXX = g++

all:$(out)

export_beta: export_beta.cpp angle_info.cpp angle_info.h ../include/sm150_parameter_paths.h
	$(GXX) $(CFLAGS) $(ROOTCFLAGS) -o $@ export_beta.cpp angle_info.cpp $(ROOTLIBS)

.PHONY: all test clean

test: tests/beta_parameters
	./tests/beta_parameters .. ../include/parameter/doppler_beta_provenance.json ../include/parameter/doppler.json

tests/beta_parameters: tests/beta_parameters.cpp angle_info.cpp angle_info.h ../include/sm150_parameter_paths.h
	$(GXX) $(CFLAGS) $(ROOTCFLAGS) -I. -o $@ tests/beta_parameters.cpp angle_info.cpp $(ROOTLIBS)

$(out): $(sourcefile)
	$(GXX) $(CFLAGS) $(ROOTCFLAGS) -o $@ $(filter %.cpp ,$(sourcefile)) $(ROOTLIBS)

clean:
	rm -f *~ *.o $(out)
