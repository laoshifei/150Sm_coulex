#ifndef SM150_GE_WALK_CANDIDATE_H
#define SM150_GE_WALK_CANDIDATE_H
// 导出完整one-pass候选；不安装、不修改active，不重复加failed旧模型。
#include "sm150_event_calibration.h"
#include "sm150_time_alignment.h"
#include <functional>
#include <set>
namespace sm150 { namespace diagnostic {
using ModelEvaluator=std::function<double(const json&,double)>;
inline json merge_log_curve(const json& old,const json& model,const ModelEvaluator& evaluate) {
  auto x=old.at("knots_keV").get<std::vector<double>>(),y=old.at("offset_ns").get<std::vector<double>>();
  require(x.size()==y.size() && x.size()>=2 && x.front()==40. && x.back()==2000.,"Invalid baseline curve");
  auto baseline=[&](double e) {
    if(e<=x.front())return y.front();
    if(e>=x.back())return y.back();
    size_t i=std::upper_bound(x.begin(),x.end(),e)-x.begin()-1;
    double f=std::log(e/x[i])/std::log(x[i+1]/x[i]);return y[i]*(1-f)+y[i+1]*f;
  };
  std::set<double> starts(x.begin(),x.end());
  for(double e:model.at("support_keV").get<std::vector<double>>())starts.insert(e);
  for(double e:model.value("knots_keV",std::vector<double>{}))starts.insert(e);
  std::vector<double> outx,outy;
  auto total=[&](double e){return baseline(e)+evaluate(model,e);};
  std::function<void(double,double,int)> interval=[&](double a,double b,int depth) {
    double ya=total(a),yb=total(b),worst=0;
    for(double f:{.25,.5,.75}) {double e=std::exp((1-f)*std::log(a)+f*std::log(b));worst=std::max(worst,std::abs(total(e)-((1-f)*ya+f*yb)));}
    if(worst>.001) {require(depth<16,"Candidate interpolation did not converge");double m=std::sqrt(a*b);interval(a,m,depth+1);interval(m,b,depth+1);}
    else {outx.push_back(b);outy.push_back(yb);}
  };
  outx.push_back(*starts.begin());outy.push_back(total(*starts.begin()));
  for(auto a=starts.begin(),b=std::next(a);b!=starts.end();++a,++b)interval(*a,*b,0);
  json merged=old;merged["knots_keV"]=outx;merged["offset_ns"]=outy;
  merged["candidate_residual_model"]=model;merged["candidate_interpolation_tolerance_ns"]=.001;
  merged["candidate_status"]="internal_validation_only_not_adopted";return merged;
}
inline json export_complete_candidate(const json& result,const json& input_meta,const fs::path& output,const ModelEvaluator& evaluate) {
  std::vector<json> branches;for(const auto& b:result.at("branches"))if(b.at("candidate_exportable").get<bool>())branches.push_back(b);
  if(branches.empty())return {{"status","not_created_no_supported_candidate"}};
  std::set<std::string> runs;
  for(const auto& input:result.at("inputs")) {
    if(!input.at("source").at("full_run").get<bool>())return {{"status","not_created_partial_input"},{"note","prefixes are code checks, not full-run parameter evidence"}};
    runs.insert(input.at("run").get<std::string>());
  }
  require(!runs.empty(),"Candidate lacks source runs");
  // 当前生产基线由提取时保存的完整路径和SHA确定，不按active重新猜测。
  fs::path source;
  for(auto it=input_meta.at("parameter_snapshots").begin();it!=input_meta.at("parameter_snapshots").end();++it)
    if(fs::path(it.key()).filename()=="calibration.json" && it.value().at("sha256")==input_meta.at("profile_sha256"))source=it.key();
  require(!source.empty(),"Cannot locate frozen complete profile");
  auto original_text=read_text(source);require(hash(original_text)==input_meta.at("profile_sha256"),"Frozen profile changed before candidate export");
  auto data=json::parse(original_text);const auto original=data;
  json changed=json::array();std::map<std::string,json> models;
  for(const auto& r:runs)for(const auto& b:branches) {
    const int ge=b.at("ge"),state=b.at("state")=="failed";auto model=b.at("selection").at("selected");
    require(!model.is_null(),"Exportable branch has no model");
    const std::string key=r+":"+std::to_string(ge)+":"+std::to_string(state);
    data["ge_curves"][key]=merge_log_curve(data.at("ge_curves").at(key),model,evaluate);models[key]=model;
    changed.push_back({{"run",r},{"ge",ge},{"state",state},{"primary_reference",b.at("primary_reference")},{"model",model}});
  }
  const fs::path directory=output/"candidate_profile";require(!fs::exists(directory),"Candidate directory exists");fs::create_directories(directory);
  data["version"]="ge_walk_candidate_"+result.at("group").get<std::string>();
  data["candidate_export"]={{"production_active",false},{"synthetic_test_only",result.value("synthetic_test_only",false)},
    {"baseline_sha256",input_meta.at("profile_sha256")},{"updated",changed},{"run_scope","only runs fully covered by this fit; other entries inherited"},
    {"formula","complete new F = complete old F + residual; existing failed_center_models retained exactly once"}};
  for(const auto& key:{"particle_tables","particle_calibration","channel_offsets","failed_center_models"})require(data.at(key)==original.at(key),"Unintended non-Ge change");
  std::string contents=data.dump();std::string identity=hash(contents);
  auto write=[&](const fs::path& p,const std::string& text){std::ofstream f(p);require(bool(f),"Cannot create candidate file");f<<text;f.close();require(bool(f),"Cannot write candidate file");};
  write(directory/"calibration.json",contents);write(directory/"source_sha256.txt",identity+"\n");
  write(directory/"correction_mode.txt","total\n");write(directory/"storage_format.txt","json\n");
  auto metadata=json::parse(read_text(source.parent_path()/"metadata.json"));metadata["candidate_export"]=data["candidate_export"];
  metadata["version"]=data["version"];write(directory/"metadata.json",metadata.dump(2)+"\n");
  double worst=0;uint64_t checked=0;
  // 用实际生产C++读取器复算，不只检查导出的JSON格式。
  for(const auto& r:runs) {
    sm150::TimeAlignment old_t,new_t;old_t.LoadProfile(source.parent_path().string(),std::stoi(r));new_t.LoadProfile(directory.string(),std::stoi(r));
    for(int ge:data.at("ge").get<std::vector<int>>())for(int state=0;state<2;++state)for(int i=-1;i<=4001;++i) {
      double e=i<0?39.:i>4000?2000.:40.*std::exp(i*std::log(1999.999/40.)/4000.);
      std::string key=r+":"+std::to_string(ge)+":"+std::to_string(state);double delta=0;
      if(e>=40 && e<2000 && models.count(key))delta=evaluate(models.at(key),e);
      double error=std::abs(new_t.TotalGeOffset(ge,state,e)-old_t.TotalGeOffset(ge,state,e)-delta);
      worst=std::max(worst,error);require(error<.002,"One-pass candidate readback mismatch");++checked;
    }
  }
  return {{"status","complete_candidate_not_adopted"},{"directory",fs::absolute(directory).string()},{"sha256",identity},{"updated_branches",changed},
    {"readback_points",checked},{"maximum_readback_error_ns",worst},{"production_active",false},{"synthetic_test_only",result.value("synthetic_test_only",false)}};
}
}}
#endif
