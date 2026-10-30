#ifndef SM150_S3_COIN_ANALYSIS_H
#define SM150_S3_COIN_ANALYSIS_H
#include <cstdint>
#include <filesystem>
#include <string>
#include "set.h"

struct AnalysisOptions {
  int run=0, segments=20;
  int64_t sample_entries=0;
  std::filesystem::path input, output, project, energy_file;
  std::string format="auto";
  double energy_min=2500., energy_max=100000.;
  double face_time=80., face_energy=2500., particle_time=80.;
  double phi_low=160., phi_high=200., cm_low=160., cm_high=200.;
  double beam=EBeam, deficit=2., projectile_mass=APro, target_mass=ATar;
};

class analysis {
 public:
  explicit analysis(AnalysisOptions options):options_(std::move(options)){}
  int Process();
 private:
  AnalysisOptions options_;
};
#endif
