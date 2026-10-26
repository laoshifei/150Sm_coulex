#ifndef SM150_GE_WALK_BUILD_H
#define SM150_GE_WALK_BUILD_H
#include <cstdint>
#include <filesystem>
#include <string>
// 标定诊断，不写活动参数；时间方向统一为particle-Ge。
struct WalkOptions {
  int run=0;
  std::filesystem::path project,input,output,profile_dir;
  std::string profile;
  uint64_t begin=0,entries=0,chunk=100000,max_buffer_hits=2000000;
  double native_window_ns=2000.,time_range_ns=1000.,time_bin_ns=4.,energy_bin_keV=10.;
  double ge_min=40.,ge_max=2000.,si_min=15000.,si_max=100000.,s3_min=25000.,s3_max=100000.;
  double s3_match_ns=80.,s3_match_keV=2500.;
  bool write_pairs=false,inspect_parameters=false;
};
void build_ge_time_walk(const WalkOptions& options);
#endif
