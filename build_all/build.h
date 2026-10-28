#ifndef SM150_BUILD_H
#define SM150_BUILD_H

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>

// Runtime settings for one build run. Empty paths use the project defaults in build.cpp.
// Read energy calibration, timing corrections and per-Ge coincidence windows from existing parameter files.
struct BuildOptions {
  int run = 0;
  std::filesystem::path project, input, outdir, receipt;
  std::filesystem::path window_file, profile_dir, source_manifest, cache_receipt;
  // windows: current narrow/wide sets; cache: wide-window raw fields; replay: rebuild narrow/wide sets from the cache.
  std::string mode = "windows", tag;

  // 0 processes the full run; a nonzero value limits a diagnostic sample, with coverage recorded in the receipt.
  uint64_t entry_limit = 0;
  std::size_t chunk_entries = 1000000;
  bool inspect_parameters = false;
};

class build {
public:
  explicit build(const BuildOptions &options);

  // Read -> energy/timing corrections -> grouping/matching -> narrow/wide gates -> output -> full readback.
  // Return 0 when processing and readback validation finish; return 1 on failure and retain error records.
  int Process();

private:
  BuildOptions options_;
};

#endif
