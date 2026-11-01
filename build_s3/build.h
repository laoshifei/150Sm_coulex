#ifndef SM150_BUILD_S3_H
#define SM150_BUILD_S3_H
#include <cstdint>
#include <filesystem>
// Unmerged calibration events before front/back normalization.
struct BuildOptions {
  int run=0;
  std::filesystem::path project,input,output,profile_dir;
  int64_t entries=0,cluster_gap_ns=100,reorder_ns=20000;
  int64_t sample_entries=-1;
  double min_energy_keV=2500.;
  int blocks=40;
};
class build {
 public:
  explicit build(BuildOptions options):options_(options){}
  int Process();
  int Sample();
 private:
  BuildOptions options_;
};
#endif
