#include "../include/sm150_source_calibration.h"

int compute_average(const sm150::sourcecal::Options& options) {
  using namespace sm150::sourcecal;
  require(!options.keys.empty(),"Average estimates require explicit --keys; there is no automatic replacement");
  const auto input=json::parse(read_text(options.input));
  require(input.at("schema")=="sm150_source_energy_candidate_v2","Expected source calibration JSON");
  json estimates=json::object();
  for(int key:options.keys) {
    const int sid=key%10000/100;require(key/10000==1,"Expected particle electronics keys");
    require(input.at("channels").contains(std::to_string(key)),"Unknown estimate target");
    std::vector<Coeff> values;std::vector<int> donors;
    for(auto it=input.at("channels").begin();it!=input.at("channels").end();++it) {
      const int other=std::stoi(it.key());
      if(other/10000!=1||other%10000/100!=sid||options.keys.count(other)||!it.value().value("candidate_eligible",false))continue;
      donors.push_back(other);values.push_back(it.value().at("calibration").at("coefficients").get<Coeff>());
    }
    json result={{"key",key},{"donor_keys",donors},{"status","insufficient_donors"},{"independent_calibration",false}};
    if(values.size()>=3) {
      Coeff mean{},spread{};for(const auto& v:values)for(int i=0;i<3;++i)mean[i]+=v[i]/values.size();
      for(const auto& v:values)for(int i=0;i<3;++i)spread[i]+=std::pow(v[i]-mean[i],2)/(values.size()-1);
      for(double& x:spread)x=std::sqrt(x);
      result["status"]="same_board_estimate_not_installed";result["coefficients"]=mean;
      result["donor_sample_stddev"]=spread;result["fit_chi2"]=nullptr;
      result["uncertainty_scope"]="Donor scatter is descriptive, not a measured uncertainty for the target channel";
    }
    estimates[std::to_string(key)]=result;
  }
  fs::create_directories(options.output);
  write_json(options.output/"estimates.json",{{"schema","sm150_same_board_estimates_v1"},
    {"input",options.input.string()},{"input_sha256",sha_file(options.input)},{"formal_replacement",false},{"estimates",estimates}});
  return 0;
}

int remap_low(const sm150::sourcecal::Options& options) {
  using namespace sm150::sourcecal;
  auto input=json::parse(read_text(options.input));
  require(input.at("schema")=="sm150_source_energy_candidate_v2"&&input.at("run")=="0457","Only the documented 0457 source remapping is implemented");
  json rows=json::object(),mapping=json::array();
  for(auto it=input.at("channels").begin();it!=input.at("channels").end();++it) {
    int old=std::stoi(it.key()),sid=old%10000/100,next=sid;
    if(sid==2)next=9;else if(sid==9)next=2;else if(sid==3)next=10;else if(sid==10)next=3;
    const int key=10000+100*next+old%100;require(!rows.contains(std::to_string(key)),"Remapping collision");
    auto row=it.value();row["key"]=key;row["source_key"]=old;rows[std::to_string(key)]=row;
    mapping.push_back({{"source_key",old},{"low_beam_key",key}});
  }
  input["channels"]=rows;input["remapping"]={{"source",options.input.string()},{"sha256",sha_file(options.input)},
    {"rule","sid2<->9 and sid3<->10; coefficient rows only; hardware gains remain separate"},{"mapping",mapping}};
  input["derived_low_beam_layout"]=true;fs::create_directories(options.output);write_json(options.output/"calibration.json",input);
  return 0;
}
