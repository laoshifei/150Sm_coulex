// 逐Ge/CFD的能量分段核心拟合；C++计算，所有未通过项保留。
// 输入为新build产生的整数二维谱，不使用ProfileY，不自动安装候选。
#include "sm150_event_calibration.h"
#include "candidate.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TObjString.h>
#include <TMatrixD.h>
#include <TDecompSVD.h>
#include <TVectorD.h>
#include <TROOT.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
using namespace sm150::diagnostic;
namespace {
struct Options {fs::path inputs,output;std::string group;int ge=0;double halfwidth=80,min_counts=100;};
const std::vector<double> ENERGY={40,60,80,100,120,160,200,250,300,370,440,520,650,800,1000,1250,1600,2000};
const std::vector<int> GE={204,205,206,209,211,212,213,214,215,304,305,306,308,310,311,506,507,508,509,510};
json fit_core(TH1D& h,double halfwidth,double min_counts) {
  json r={{"status","insufficient_counts"},{"mu_ns",nullptr},{"error_ns",nullptr},{"counts",h.Integral()}};
  if(h.Integral()<min_counts) return r;
  int peak=0;double best=-1;const double bw=h.GetBinWidth(1);
  for(int i=3;i<=h.GetNbinsX()-2;++i) if(std::abs(h.GetBinCenter(i))<250) {
    double v=h.GetBinContent(i-2)+2*h.GetBinContent(i-1)+3*h.GetBinContent(i)+2*h.GetBinContent(i+1)+h.GetBinContent(i+2);
    if(v>best) {peak=i;best=v;}
  }
  if(!peak) {r["status"]="no_search_bins";return r;}
  const double seed=h.GetBinCenter(peak);
  if(std::abs(seed)>238) {r["status"]="peak_at_search_edge";r["seed_ns"]=seed;return r;}
  std::vector<double> bg;
  for(int i=1;i<=h.GetNbinsX();++i) if(std::abs(h.GetBinCenter(i)-seed)>2*halfwidth) bg.push_back(h.GetBinContent(i));
  std::sort(bg.begin(),bg.end());double background=bg.empty()?0:bg[bg.size()/2];
  double height=h.GetBinContent(peak);
  if(height-background<3*std::sqrt(std::max(1.,background))) {r["status"]="no_clear_peak";return r;}
  const double lo=seed-halfwidth,hi=seed+halfwidth;
  TF1 fun((std::string(h.GetName())+"_core").c_str(),[lo,hi](double* x,double* p){double u=(x[0]-p[1])/p[2];return p[0]*std::exp(-.5*u*u)+p[3]+(p[4]-p[3])*(x[0]-lo)/(hi-lo);},lo,hi,5);
  fun.SetParameters(std::max(1.,height-background),seed,16.,background,background);
  fun.SetParLimits(0,0,std::max(20.,height*20));fun.SetParLimits(1,lo+2*bw,hi-2*bw);
  fun.SetParLimits(2,4.,std::min(100.,halfwidth/1.5));
  fun.SetParLimits(3,0,std::max(20.,height*5));fun.SetParLimits(4,0,std::max(20.,height*5));
  TFitResultPtr fit=h.Fit(&fun,"QRS0LNI");int status=int(fit),cov=fit.Get()?fit->CovMatrixStatus():-1;
  double mu=fun.GetParameter(1),error=fun.GetParError(1),sigma=fun.GetParameter(2);
  r["root_status"]=status;r["covariance_status"]=cov;r["fit_range_ns"]={lo,hi};r["sigma_ns"]=sigma;
  r["root_fit_statistic"]=fun.GetChisquare();r["ndf"]=fun.GetNDF();r["parameters"]={fun.GetParameter(0),mu,sigma,fun.GetParameter(3),fun.GetParameter(4)};
  r["fitted_mu_ns"]=std::isfinite(mu)?json(mu):json(nullptr);r["fitted_error_ns"]=std::isfinite(error)?json(error):json(nullptr);
  r["core_area"]=fun.GetParameter(0)*sigma*std::sqrt(2*std::acos(-1.))/bw;
  if(status!=0 || cov!=3 || !std::isfinite(error) || error<=0 || error>20 || sigma<=4.01 || sigma>=std::min(100.,halfwidth/1.5)*.995 || mu<=lo+2.01*bw || mu>=hi-2.01*bw) {r["status"]="fit_failed_or_boundary";return r;}
  if(fun.GetNDF()<=0 || fun.GetChisquare()/fun.GetNDF()>5) {r["status"]="core_shape_mismatch";return r;}
  if(r["core_area"].get<double>()<25 || fun.GetParError(0)<=0 || fun.GetParameter(0)/fun.GetParError(0)<4) {r["status"]="weak_core";return r;}
  for(int i=3;i<=h.GetNbinsX()-2;++i) if(std::abs(h.GetBinCenter(i))<250 && std::abs(h.GetBinCenter(i)-mu)>halfwidth) {
    double v=h.GetBinContent(i-2)+2*h.GetBinContent(i-1)+3*h.GetBinContent(i)+2*h.GetBinContent(i+1)+h.GetBinContent(i+2);
    if(v-9*background>.8*(best-9*background)) {r["status"]="ambiguous_secondary_peak";return r;}
  }
  r["status"]="valid";r["mu_ns"]=mu;r["error_ns"]=error;return r;
}
json nodes(TH2D& h,const Options& o,const std::string& label) {
  json answer=json::array();
  for(size_t i=0;i+1<ENERGY.size();++i) {
    double lo=ENERGY[i],hi=ENERGY[i+1];int first=h.GetXaxis()->FindBin(lo+1e-6),last=h.GetXaxis()->FindBin(hi-1e-6);
    std::unique_ptr<TH1D> slice(h.ProjectionY((label+"_"+std::to_string(i)).c_str(),first,last));slice->SetDirectory(nullptr);
    auto r=fit_core(*slice,o.halfwidth,o.min_counts);r["energy_bin"]=i;r["energy_lo_keV"]=lo;r["energy_hi_keV"]=hi;r["energy_keV"]=(lo+hi)/2;
    if(r["status"]=="valid") {
      double mu=r["mu_ns"],sigma=r["sigma_ns"],sum=0,sumE=0;
      int y0=h.GetYaxis()->FindBin(mu-2*sigma),y1=h.GetYaxis()->FindBin(mu+2*sigma);
      for(int x=first;x<=last;++x) for(int y=y0;y<=y1;++y) {double n=h.GetBinContent(x,y);sum+=n;sumE+=n*h.GetXaxis()->GetBinCenter(x);}
      if(sum>0) r["energy_keV"]=sumE/sum;
      auto check=fit_core(*slice,o.halfwidth*1.25,o.min_counts);r["window_check"]=check;
      if(check["status"]!="valid") r["status"]="window_check_unresolved";
      else if(std::abs(check["mu_ns"].get<double>()-mu)>std::max(4.,3*std::hypot(check["error_ns"].get<double>(),r["error_ns"].get<double>()))) r["status"]="window_sensitive";
    }
    answer.push_back(r);
  }
  return answer;
}
std::vector<double> basis(const std::string& family,double e,const std::vector<double>& knots) {
  double u=std::log(e/334.);
  if(family=="constant") return {1};
  if(family=="inverse_sqrt") return {1,std::sqrt(334./e)-1};
  if(family=="inverse") return {1,334./e-1};
  if(family=="log_linear") return {1,u};
  if(family=="log_quadratic") return {1,u,u*u};
  std::vector<double> b(knots.size(),0.);
  if(e<=knots.front()) b.front()=1.;
  else if(e>=knots.back()) b.back()=1.;
  else {size_t i=std::upper_bound(knots.begin(),knots.end(),e)-knots.begin()-1;double f=std::log(e/knots[i])/std::log(knots[i+1]/knots[i]);b[i]=1-f;b[i+1]=f;}
  return b;
}
double value(const json& model,double e) {
  auto domain=model.at("support_keV").get<std::vector<double>>();e=std::max(domain[0],std::min(domain[1],e));
  auto b=basis(model.at("family"),e,model.value("knots_keV",std::vector<double>{}));auto c=model.at("coefficients_ns").get<std::vector<double>>();
  double y=0;for(size_t i=0;i<b.size();++i)y+=b[i]*c[i];return y;
}
json curve(const json& points,const std::string& family) {
  std::vector<double> e,y,error;
  for(const auto& p:points) if(p.at("status")=="valid") {e.push_back(p.at("energy_keV"));y.push_back(p.at("mu_ns"));error.push_back(std::max(1.,p.at("error_ns").get<double>()));}
  if(e.empty()) return nullptr;
  std::vector<double> knots;
  if(family=="piecewise_log_linear") for(int i=0;i<4;++i) knots.push_back(std::exp(std::log(e.front())+i*(std::log(e.back())-std::log(e.front()))/3));
  const int k=basis(family,e.front(),knots).size(),n=e.size();if(n<k+2 || (family!="constant" && e.back()/e.front()<1.8)) return nullptr;
  TMatrixD a(n,k);TVectorD b(n);
  for(int i=0;i<n;++i) {auto row=basis(family,e[i],knots);b[i]=y[i]/error[i];for(int j=0;j<k;++j)a(i,j)=row[j]/error[i];}
  TDecompSVD svd(a);Bool_t ok;auto coeff=svd.Solve(b,ok);if(!ok) return nullptr;
  auto singular=svd.GetSig();double condition=singular[0]/singular[k-1];if(!std::isfinite(condition)||condition>1e8)return nullptr;
  std::vector<double> c(k);for(int i=0;i<k;++i)c[i]=coeff[i];
  json m={{"family",family},{"coefficients_ns",c},{"support_keV",{e.front(),e.back()}},{"knots_keV",knots},{"n_points",n},{"n_parameters",k},{"svd_condition",condition},{"weight_floor_ns",1.},{"endpoint_policy","clamp_new_residual_to_measured_support"}};
  double ss=0;for(int i=0;i<n;++i)ss+=std::pow(y[i]-value(m,e[i]),2);m["fit_rms_ns"]=std::sqrt(ss/n);return m;
}
json evaluate(const json& model,const json& check) {
  double before=0,after=0;int n=0;auto domain=model.at("support_keV").get<std::vector<double>>();
  for(const auto& p:check) if(p.at("status")=="valid") {
    double e=p.at("energy_keV");if(e<domain[0] || e>domain[1])continue;
    double y=p.at("mu_ns");before+=y*y;after+=std::pow(y-value(model,e),2);++n;
  }
  if(n<3) return nullptr;
  return {{"n",n},{"before_rms_ns",std::sqrt(before/n)},{"after_rms_ns",std::sqrt(after/n)}};
}
int supported(const json& points) {
  int n=0;for(const auto& p:points) if(p.at("status")=="valid") ++n;return n;
}
// 两个时间半段互相检验，不把同一事件的Si/S3参考称为独立观测。
json select_model(const json& full,const json& first,const json& second) {
  json answer={{"status","insufficient_support"},{"selected",nullptr},{"models",json::array()}};
  double best=std::numeric_limits<double>::infinity();int best_k=100;
  for(const std::string family:{"constant","inverse_sqrt","inverse","log_linear","log_quadratic","piecewise_log_linear"}) {
    auto m=curve(full,family),a=curve(first,family),b=curve(second,family);
    if(m.is_null()) continue;
    json item={{"full_model",m},{"first_model",a},{"second_model",b},{"passes_halves",false}};
    if(!a.is_null() && !b.is_null()) {
      auto ab=evaluate(a,second),ba=evaluate(b,first);item["first_to_second"]=ab;item["second_to_first"]=ba;
      // 合并谱可能失去部分能量点，最终曲线不一定等于两半各自的好曲线。
      // 必须检查实际导出的完整模型，防止在缺点能区出现失控插值。
      auto ma=evaluate(m,first),mb=evaluate(m,second);
      item["full_model_on_first"]=ma;item["full_model_on_second"]=mb;
      if(!ab.is_null() && !ba.is_null() && !ma.is_null() && !mb.is_null()) {
        const double score=std::sqrt((std::pow(ab["after_rms_ns"].get<double>(),2)+std::pow(ba["after_rms_ns"].get<double>(),2))/2);
        bool pass=true;
        for(const auto& ev:{ab,ba,ma,mb}) {
          const double old=ev.at("before_rms_ns"),now=ev.at("after_rms_ns");
          if(old-now<std::max(1.,.1*old)) pass=false;
        }
        item["passes_halves"]=pass;item["validation_rms_ns"]=score;
        const int k=m.at("n_parameters");
        if(pass && (score<best-1. || (std::abs(score-best)<=1. && k<best_k))) {
          answer["selected"]=m;answer["validation_rms_ns"]=score;answer["status"]="candidate_halves_passed";best=score;best_k=k;
        }
      }
    }
    answer["models"].push_back(item);
  }
  if(answer["selected"].is_null() && !answer["models"].empty()) answer["status"]="keep_current_no_repeated_improvement";
  return answer;
}
json other_reference(const json& model,const json& points) {
  std::vector<double> old,now;auto support=model.at("support_keV").get<std::vector<double>>();
  for(const auto& p:points) if(p.at("status")=="valid") {
    double e=p.at("energy_keV");if(e<support[0] || e>support[1]) continue;
    old.push_back(p.at("mu_ns"));now.push_back(old.back()-value(model,e));
  }
  if(old.size()<3) return {{"status","insufficient_other_reference"},{"n",old.size()}};
  double mean0=0,mean1=0;for(size_t i=0;i<old.size();++i){mean0+=old[i];mean1+=now[i];}mean0/=old.size();mean1/=old.size();
  double ss0=0,ss1=0;for(size_t i=0;i<old.size();++i){ss0+=std::pow(old[i]-mean0,2);ss1+=std::pow(now[i]-mean1,2);}
  const double rms0=std::sqrt(ss0/old.size()),rms1=std::sqrt(ss1/old.size());
  return {{"status",rms1>rms0+2.?"shape_disagreement":"shape_not_worse"},{"n",old.size()},
    {"before_mean_ns",mean0},{"after_mean_ns",mean1},{"before_centered_rms_ns",rms0},{"after_centered_rms_ns",rms1},
    {"note","constant Si-S3 reference difference is reported, not absorbed into a Ge energy term"}};
}
void save(const fs::path& path,const json& data) {
  auto temporary=path.string()+".tmp";std::ofstream f(temporary);require(bool(f),"Cannot create "+temporary);
  f<<data.dump(2)<<'\n';f.close();require(bool(f),"Cannot write "+temporary);fs::rename(temporary,path);
}
std::string name(int ge,const std::string& ref,const std::string& state,int half,const std::string& coordinate) {
  return "ge"+std::to_string(ge)+"_"+ref+"_"+state+"_half"+std::to_string(half)+"_"+coordinate;
}
void run_fit(const Options& o) {
  require(o.group=="low" || o.group=="high1" || o.group=="high2","Require --group low|high1|high2");
  require(o.ge==0 || std::find(GE.begin(),GE.end(),o.ge)!=GE.end(),"Unsupported Ge");
  require(std::isfinite(o.halfwidth) && std::isfinite(o.min_counts) && o.halfwidth>=40 && o.halfwidth<=120 && o.min_counts>=25,"Invalid fit settings");
  require(!fs::exists(o.output),"Output directory exists; preserve earlier fits");
  auto paths=json::parse(read_text(o.inputs));if(paths.is_object()) paths=paths.at("inputs");
  require(paths.is_array() && !paths.empty(),"Inputs must be a nonempty path array");
  std::map<std::string,std::unique_ptr<TH2D>> summed;
  std::map<std::string,std::vector<std::pair<uint64_t,uint64_t>>> intervals;
  json sources=json::array(),baseline,signature;std::set<std::string> seen;
  for(const auto& path:paths) {
    fs::path p=fs::canonical(path.get<std::string>());require(seen.insert(p.string()).second,"Duplicate ROOT input");
    Stat before(p);std::unique_ptr<TFile> f(TFile::Open(p.c_str(),"READ"));require(f && !f->IsZombie(),"Cannot read "+p.string());
    auto* text=dynamic_cast<TObjString*>(f->Get("metadata_json"));require(text,"Missing metadata");auto m=json::parse(text->GetString().Data());
    require(m.at("schema")=="sm150_ge_walk_v1" && m.at("status")=="complete","Incomplete or incompatible ROOT");
    require(m.at("group")==o.group,"Mixed periods are not allowed");
    require(m.at("energy_bin_keV")==10. && m.at("energy_gates_keV").at("ge")==json::array({40.,2000.}),"Fit uses fixed40-2000keV,10keV source bins");
    json sig;for(const auto& key:{"schema","timing_profile","profile_sha256","energy_gates_keV","native_pair_window_ns","time_range_ns","time_bin_ns","energy_bin_keV","dt_sign","s3_match_ns","s3_match_keV","quality"})sig[key]=m.at(key);
    if(signature.is_null()) {signature=sig;baseline=m;} else require(signature==sig,"Mismatched extraction configuration");
    const std::string run=m.at("run");uint64_t begin=m.at("source").at("entry_begin"),end=m.at("source").at("entry_end");
    for(auto old:intervals[run])require(end<=old.first || begin>=old.second,"Overlapping input intervals");
    intervals[run].push_back({begin,end});
    for(int g:GE) if(!o.ge || g==o.ge) for(const std::string ref:{"si","s3_matched"})
      for(const std::string state:{"normal","failed"}) for(int half=0;half<2;++half) for(const std::string coordinate:{"current","native_ge"}) {
        auto n=name(g,ref,state,half,coordinate);auto* h=dynamic_cast<TH2D*>(f->Get(n.c_str()));require(h,"Missing histogram "+n);
        require(h->GetNbinsX()==196 && h->GetXaxis()->GetXmin()==40. && h->GetXaxis()->GetXmax()==2000.,"Unexpected energy edges");
        require(h->GetNbinsY()*m.at("time_bin_ns").get<double>()==2*m.at("time_range_ns").get<double>(),"Unexpected time bins");
        if(!summed.count(n)) {auto* c=dynamic_cast<TH2D*>(h->Clone());c->SetDirectory(nullptr);summed[n].reset(c);}
        else require(summed[n]->Add(h),"Incompatible histograms");
    }
    require(before.same(Stat(p)),"Input changed during fitting preparation");
    sources.push_back({{"path",p.string()},{"sha256",hash(read_text(p))},{"run",run},{"source",m.at("source")},{"extractor_sha256",m.at("program_sha256")}});
  }
  fs::create_directories(o.output);
  json result={{"schema","sm150_ge_walk_core_fit_v1"},{"status","running"},{"group",o.group},{"inputs",sources},{"extraction",signature},
    {"fitter_sha256",hash(read_text(fs::canonical("/proc/self/exe")))},{"ROOT_version",gROOT->GetVersion()},
    {"method",{{"slice_model","Gaussian core + nonnegative linear background, Poisson bin-integrated fit"},{"core_halfwidth_ns",o.halfwidth},{"window_check_factor",1.25},{"minimum_slice_counts",o.min_counts},{"energy_edges_keV",ENERGY},
      {"energy_coordinate","source10keV-bin centers weighted by observed counts within fitted core mu+-2sigma"},
      {"validation","first/second temporal halves, both directions; internal consistency, not blind validation"},
      {"selection","both cross-half directions and the final full-data curve on each half improve RMS by max(1ns,10%); prefer simpler within1ns CV RMS"},
      {"fit_error_limit","counts from shared hits can be correlated; formal Poisson errors are conditional diagnostics"},
      {"parameter_sign","new Ge full correction = existing Ge full correction + fitted current(particle-Ge) residual"},
      {"production_active",false}}},{"branches",json::array()}};
  const auto start=std::chrono::steady_clock::now();
  for(int g:GE) if(!o.ge || g==o.ge) for(const std::string state:{"normal","failed"}) {
    json branch={{"ge",g},{"state",state},{"references",json::object()}};
    for(const std::string ref:{"s3_matched","si"}) {
      auto& a=*summed.at(name(g,ref,state,0,"current"));auto& b=*summed.at(name(g,ref,state,1,"current"));
      std::string stem="ge"+std::to_string(g)+"_"+ref+"_"+state;
      std::unique_ptr<TH2D> full(dynamic_cast<TH2D*>(a.Clone((stem+"_full").c_str())));full->SetDirectory(nullptr);full->Add(&b);
      json n={{"full",nodes(*full,o,stem+"_full")},{"first",nodes(a,o,stem+"_first")},{"second",nodes(b,o,stem+"_second")},{"entries",full->GetEntries()}};
      n["selection"]=select_model(n["full"],n["first"],n["second"]);branch["references"][ref]=n;
    }
    auto& refs=branch["references"];
    std::string primary="s3_matched";
    if(supported(refs[primary]["full"])<4 || supported(refs[primary]["first"])<3 || supported(refs[primary]["second"])<3) primary="si";
    const std::string secondary=primary=="si"?"s3_matched":"si";
    branch["primary_reference"]=primary;branch["selection"]=refs[primary]["selection"];
    auto selected=branch["selection"]["selected"];branch["candidate_exportable"]=false;
    if(!selected.is_null()) {
      auto cross=other_reference(selected,refs[secondary]["full"]);branch["other_reference_check"]=cross;
      if(cross["status"]=="shape_disagreement") branch["selection"]["status"]="reference_shape_disagreement";
      else {branch["candidate_exportable"]=true;if(cross["status"]=="insufficient_other_reference")branch["selection"]["status"]="candidate_halves_passed_single_reference";}
    }
    result["branches"].push_back(branch);save(o.output/"fit_results.json",result);
    std::cout<<g<<' '<<state<<" reference="<<primary<<" nodes="<<supported(refs[primary]["full"])<<" "<<branch["selection"]["status"]<<'\n'<<std::flush;
  }
  {TFile out((o.output/"aggregate.root").c_str(),"CREATE");require(!out.IsZombie(),"Cannot create aggregate ROOT");
   for(const auto& item:summed)require(item.second->Write()>0,"Cannot save histogram");
   require(!out.TestBit(TFile::kWriteError),"Aggregate ROOT write failed");out.Close();}
  std::ofstream table(o.output/"branches.csv");table<<"group,ge,state,reference,status,valid_nodes,candidate_exportable,family,support_lo_keV,support_hi_keV\n";
  for(const auto& b:result["branches"]) {
    std::string ref=b.at("primary_reference");auto m=b.at("selection").at("selected");
    table<<o.group<<','<<b.at("ge")<<','<<b.at("state").get<std::string>()<<','<<ref<<','<<b.at("selection").at("status").get<std::string>()<<','<<supported(b.at("references").at(ref).at("full"))<<','<<b.at("candidate_exportable")<<',';
    if(!m.is_null())table<<m.at("family").get<std::string>()<<','<<m.at("support_keV")[0]<<','<<m.at("support_keV")[1];else table<<",,";table<<'\n';
  }
  table.close();require(bool(table),"Cannot write branch table");
  result["complete_candidate"]=export_complete_candidate(result,baseline,o.output,value);
  result["seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();result["status"]="complete";
  save(o.output/"fit_results.json",result);
}
}
int main(int argc,char** argv) {
  try {
    Options o;
    for(int i=1;i<argc;++i) {
      std::string key=argv[i];if(key=="--help") {std::cout<<"fit --inputs paths.json --output NEW_DIRECTORY --group low|high1|high2 [--ge KEY --halfwidth 80 --min-counts 100]\nCore fits and candidate residuals only; does not change active parameters.\n";return 0;}
      require(i+1<argc,"Missing value for "+key);std::string val=argv[++i];
      if(key=="--inputs")o.inputs=val;else if(key=="--output")o.output=val;else if(key=="--group")o.group=val;
      else if(key=="--ge")o.ge=std::stoi(val);else if(key=="--halfwidth")o.halfwidth=std::stod(val);else if(key=="--min-counts")o.min_counts=std::stod(val);else throw std::runtime_error("Unknown option "+key);
    }
    require(!o.inputs.empty() && !o.output.empty(),"Specify --inputs and --output");gROOT->SetBatch(true);TH1::AddDirectory(false);run_fit(o);return 0;
  } catch(const std::exception& e) {std::cerr<<"ge_time_walk fit: "<<e.what()<<'\n';return 1;}
}
