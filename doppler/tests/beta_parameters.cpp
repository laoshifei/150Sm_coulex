// Compare precomputed beta with its frozen source model, including every
// particle ring/sector and interleaved periods. Calibration values stay in JSON.
#include "angle_info.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

using nlohmann::json;
using namespace sm150_doppler;

void require(bool condition,const char *message) {
  if(!condition)throw std::runtime_error(message);
}
json read(const char *name) {
  std::ifstream file(name);require(bool(file),"Cannot read test parameters");
  json value;file>>value;return value;
}
void same(const Correction &a,const Correction &b) {
  require(a.energy==b.energy&&a.beta==b.beta&&a.cosine==b.cosine&&
    a.theta_p==b.theta_p&&a.phi_p==b.phi_p&&a.theta_r==b.theta_r&&a.phi_r==b.phi_r&&
    a.theta_g==b.theta_g&&a.phi_g==b.phi_g&&a.expected_particle_MeV==b.expected_particle_MeV&&
    a.extrapolated==b.extrapolated&&a.stopped==b.stopped,"Correction differs from explicit per-period reference");
}
int main(int argc,char **argv) {
  try {
    require(argc==4,"Usage: beta_parameters PROJECT PROVENANCE_JSON DIRECT_BETA_JSON");
    const std::string project=std::filesystem::absolute(argv[1]);
    const auto baseline=read(argv[2]).at("source_parameters"),adopted=read(argv[3]);
    Model current(project,adopted),old(project,baseline);
    require(current.uses_tabulated_beta()&&!old.uses_tabulated_beta(),"Wrong beta representation");
    std::size_t comparisons=0;
    for(const int start:{0,1,2}) {
      Model interleaved(project,adopted);
      for(int ref=0;ref<2;++ref)for(int ring=1;ring<=(ref ? 24:8);++ring)for(int sector=1;sector<=(ref ? 32:12);++sector) {
        Pair p;p.ge=214;p.ge_ring=3;p.ge_sector=2;p.reference=ref;p.ring=ring;p.sector=sector;
        p.raw=334.0;p.particle_MeV=50.0;
        // Repeatedly revisit identical geometry in different periods. Without
        // the period in the cache key, the first period contaminates the next.
        for(const int step:{0,1,2,0,2,1}) {
          p.group=(start+step)%3;
          const auto result=interleaved.evaluate(p),expected=old.evaluate(p);
          same(result,expected);++comparisons;
          const std::array<std::string,3> periods={"high1","low","high2"};
          const double stored=adopted.at("beta_by_period").at(periods[p.group]).at(ref ? "S3":"Si").at(std::to_string(ring)).at("beta");
          require(result.beta==stored,"Stored beta was scaled again");
        }
      }
    }
    std::size_t rejected=0;
    auto reject=[&](const json &bad) {
      try {Model model(project,bad);}
      catch(const std::exception &) {++rejected;return;}
      throw std::runtime_error("Invalid parameter configuration accepted");
    };
    auto bad=adopted;bad["schema_version"]=2;reject(bad);
    bad=adopted;bad["schema_version"]=4;reject(bad);
    bad=adopted;bad["S3"]["energy_scale"]=1.;reject(bad);
    bad=adopted;bad["energy_scale_by_period"]=json::object();reject(bad);
    bad=adopted;bad["beta_by_period"].erase("high1");reject(bad);
    bad=adopted;bad["beta_by_period"]["high1"]["S3"].erase("1");reject(bad);
    bad=adopted;bad["beta_by_period"]["high1"]["S3"]["1"]["beta"]=-1.;reject(bad);
    bad=adopted;bad["beta_by_period"]["high1"]["S3"]["1"]["beta"]=1.;reject(bad);
    bad=adopted;bad["beta_by_period"]["high1"]["S3"]["1"]["theta_p_deg"]=0.;reject(bad);
    bad=adopted;bad["S3"]["theta_offset_deg"]=1.;reject(bad);
    bad=adopted;bad["Si"]["particle_energy_fraction"]=1.;reject(bad);
    bad=adopted;bad["beta_input_sha256"]["beta_info/theta_ein_eout_r.txt"]="invalid";reject(bad);
    bool missing_scale=false;try{current.energy_scale(0,1);}catch(const std::runtime_error &){missing_scale=true;}
    require(missing_scale,"Direct beta must not expose a second scale");
    std::cout<<json({{"status","passed"},{"exact_model_comparisons",comparisons},
      {"invalid_inputs_rejected",rejected},{"geometry_unchanged",true},{"beta_applied_once",true}}).dump(2)<<'\n';
  }catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
