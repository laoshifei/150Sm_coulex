#include "../include/sm150_source_calibration.h"

int compute_k_all(const sm150::sourcecal::Options& options) {
  using namespace sm150::sourcecal;
  const auto seeds=read_seed(options.seed,0);Input input(options.input);
  require(!options.legacy,"Gain mode requires modern metadata proving energy_base has no old k");
  require(input.metadata.value("schema","")=="sm150_raw2hist_v1","Missing raw2hist metadata");
  require(input.metadata.at("run").get<std::string>()==sm150::diagnostic::rn(options.run),"Gain run mismatch");
  require(input.metadata.value("format","")=="decoded","Cannot derive total k from already-calibrated input");
  const auto digest=sha_file(options.seed);bool seed_found=false;
  for(auto it=input.metadata.at("parameter_snapshots").begin();it!=input.metadata.at("parameter_snapshots").end();++it)
    if(fs::path(it.key()).filename()=="cali_ge.txt") {
      require(it.value().at("sha256")==digest,"Histogram base calibration differs from --seed");seed_found=true;
    }
  require(seed_found,"Missing provenance for uncorrected Ge base calibration");
  fs::create_directories(options.output);
  json result={{"schema","sm150_ge_gain_candidate_v1"},{"run",sm150::diagnostic::rn(options.run)},
    {"input",options.input.string()},{"input_sha256",sha_file(options.input)},{"input_metadata",input.metadata},
    {"seed_sha256",digest},{"binary_sha256",sha_file(fs::canonical("/proc/self/exe"))},
    {"formula","k=511 keV / centroid of energy_base; old k is not applied"},
    {"formal_replacement",false},{"channels",json::object()}};
  std::ofstream csv(options.output/"gain.csv"),text(options.output/"k_candidate.txt");
  csv<<"key,status,k,k_se,old_k,center_keV,center_se_keV,deviance_per_dof\n"<<std::setprecision(17);
  text<<"cid sid ch k\n"<<std::setprecision(17);
  int selected=0,eligible=0;
  for(const auto& item:seeds) {
    const int key=item.first;if(!options.keys.empty()&&!options.keys.count(key))continue;++selected;
    json r={{"key",key},{"status","missing_histogram"},{"candidate_eligible",false}};
    if(auto* h=input.histogram(key,"energy_base",false)) {
      int at=h->FindBin(511.);
      for(int i=h->FindBin(495.);i<=h->FindBin(530.);++i)if(h->GetBinContent(i)>h->GetBinContent(at))at=i;
      r=choose_peak(*h,h->GetBinCenter(at),15.,0.);r["key"]=key;r["candidate_eligible"]=false;
      if(r.contains("center")) {
        const double center=r["center"],se=r["center_se"],k=511./center;
        const bool good=r.value("center_usable",false)&&center>495&&center<530&&
          r.at("deviance_per_dof").get<double>()<5&&k>.9&&k<1.1;
        r["k"]=k;r["k_se"]=511*se/(center*center);r["candidate_eligible"]=good;
        r["status"]=good?"gain_candidate":"gain_review_required";
        if(good){++eligible;text<<key/10000<<' '<<key%10000/100<<' '<<key%100<<' '<<k<<'\n';}
      }
    }
    const std::string name=std::to_string(key);r["old_k"]=input.metadata.at("parameter_records").at(name).at("k_ge");
    csv<<key<<','<<r.at("status").get<std::string>();
    if(r.contains("k"))csv<<','<<r["k"]<<','<<r["k_se"]<<','<<r["old_k"]<<','<<r["center"]<<','<<r["center_se"]<<','<<r["deviance_per_dof"];
    else csv<<",,,"<<r["old_k"]<<",,,";
    csv<<'\n';result["channels"][name]=r;
  }
  require(selected>0&&(options.keys.empty()||selected==int(options.keys.size())),"Unknown requested Ge keys");
  csv.close();text.close();require(bool(csv)&&bool(text),"Gain output failed");input.verify();
  result["completion"]={{"status","computation_complete"},{"selected_channels",selected},{"candidate_channels",eligible},{"candidate_text_is_complete",selected==eligible}};
  write_json(options.output/"gain.json",result);write_json(options.output/"completion.json",result["completion"]);
  return 0;
}
