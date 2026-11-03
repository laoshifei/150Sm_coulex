// Export final ring beta values through the existing correction model. No
// experimental event data are read and no numerical scale is embedded here.
#include "angle_info.h"
#include "../include/sm150_parameter_paths.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

int main(int argc,char **argv) {try {
  if(argc!=4)throw std::runtime_error("Usage: export_beta PROJECT SOURCE_JSON NEW_JSON");
  namespace fs=std::filesystem;using nlohmann::json;
  const auto project=fs::absolute(argv[1]);const fs::path source=argv[2],output=argv[3];
  if(fs::exists(output))throw std::runtime_error("Refusing to overwrite beta export");
  std::ifstream input(source);if(!input)throw std::runtime_error("Cannot read source parameters");
  json original;input>>original;
  for(const auto *family:{"Si","S3"})
    if(original.at(family).at("theta_offset_deg")!=0||original.at(family).at("particle_energy_fraction")!=0)
      throw std::runtime_error("Ring beta export requires nominal particle theta and no event-energy blend");
  sm150_doppler::Model model(project.string(),original);
  sm150::ParameterPaths paths(project/"pars");
  json result=original;
  result["schema_version"]=3;result["label"]="tabulated_beta";
  result.erase("energy_scale_by_period");
  for(const auto *family:{"Si","S3"})result[family].erase("energy_scale");
  result["beta_by_period"]=json::object();
  result["beta_units"]="v/c";
  result["source_parameters_sha256"]=sm150::ParameterPaths::sha256(source);
  result["notes"]="Final ring beta values exported from the adopted model. Apply directly without an energy or beta multiplier. Original Ge and particle angles are retained.";
  for(const auto *name:{"angle_info/spider_angle.txt","angle_info/s3_ring_angle.txt", "beta_info/theta_ein_eout_p.txt","beta_info/theta_ein_eout_r.txt"})
    result["beta_input_sha256"][name]=sm150::ParameterPaths::sha256(paths.legacy(name));
  std::size_t count=0;
  const std::array<std::string,3> periods={"high1","low","high2"};
  for(int group=0;group<3;++group)for(int ref=0;ref<2;++ref)for(int ring=1;ring<=(ref ? 24:8);++ring) {
    sm150_doppler::Pair pair;pair.group=group;pair.ge=214;pair.ge_ring=3;pair.ge_sector=2;
    pair.reference=ref;pair.ring=ring;pair.sector=1;pair.raw=334;pair.particle_MeV=50;
    const auto value=model.evaluate(pair);
    result["beta_by_period"][periods[group]][ref ? "S3":"Si"][std::to_string(ring)]={
      {"beta",value.beta},{"theta_p_deg",value.theta_p},{"theta_r_deg",value.theta_r},
      {"expected_particle_MeV",value.expected_particle_MeV},{"kinematic_extrapolated",value.extrapolated}};
    ++count;
  }
  std::ofstream file(output);if(!file)throw std::runtime_error("Cannot write beta export");
  file<<result.dump(2)<<'\n';file.close();if(!file)throw std::runtime_error("Beta export write failed");
  std::cout<<json({{"status","exported"},{"ring_instances",count},{"output",output.string()}}).dump()<<'\n';
}catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 1;}}
