#ifndef SM150_OFFSET_ANALYSIS_H
#define SM150_OFFSET_ANALYSIS_H
#include "set.h"
#include <cstdint>
#include <filesystem>
#include <string>

// Produce candidate measurements/residual diagnostics only; do not write pars or adopt offsets automatically.
struct OffsetOptions {
  int run = 0, ge_reference = KEYREF;
  std::filesystem::path project, input, output;
  std::string mode = "raw", energy_stage = "final", profile;
  double window_ns = ALIGNMENTWINDOW, bin_ns = BINWIDTH;
  double ge_min = 40., ge_max = 2000.;
  double si_min = 15000., si_max = 100000.;
  double s3_min = 25000., s3_max = 100000.;
  uint64_t begin = 0, entries = 0, chunk = 100000;
  uint64_t max_buffer_hits = 2000000;
  bool write_pairs = false, inspect_parameters = false, particle_calibration = false, particle_network = false;
  bool particle_full = false;
};
class analysis {
public:
  explicit analysis(const OffsetOptions& options) : options_(options) {}
  void Process();
private:
  OffsetOptions options_;
};
#endif
