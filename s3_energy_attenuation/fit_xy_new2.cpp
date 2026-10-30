// Relative energy-scale registration with independent spectrum normalization.
// Hellinger distances are descriptive diagnostics, not chi-square values.
#include "sm150_event_calibration.h"
#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TNamed.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <unistd.h>
using namespace sm150::diagnostic;
namespace {
constexpr int channels=56, points=241;
constexpr double bin_width=50., fit_low=20000., fit_high=90000., fit_bin=250., minimum_count=400.;
struct Spectrum {std::vector<double> fit,check;double old_k=1.;std::string name;};
struct Input {json meta;std::vector<Spectrum> spectra;std::string sha;};
Input load(const fs::path& path) {
  Input out;out.meta=json::parse(read_text(path.string()+".json"));
  require(out.meta.at("schema")=="sm150_s3_attenuation_spectra_v1" && out.meta.at("status")=="complete","Wrong source schema/status");
  out.sha=hash(read_text(path));std::unique_ptr<TFile> file(TFile::Open(path.c_str(),"READ"));
  require(file && !file->IsZombie(),"Cannot read source spectra");out.spectra.resize(channels);
  std::vector<bool> found(channels,false);
  for(const auto& c:out.meta.at("channels")) {
    int index=c.at("index");require(index>=0 && index<channels && !found[index],"Duplicate/invalid channel");found[index]=true;
    auto& s=out.spectra[index];s.name=c.at("name");s.old_k=c.at("old_k");
    auto* h=dynamic_cast<TH2D*>(file->Get(("clean_"+s.name).c_str()));require(h,"Missing spectrum "+s.name);
    require(h->GetNbinsX()==3000 && h->GetXaxis()->GetXmin()==0. && h->GetXaxis()->GetXmax()==150000.,"Wrong energy axis");
    require(h->GetNbinsY()>=4 && h->GetNbinsY()%2==0,"Wrong segment axis");
    s.fit.assign(3000,0.);s.check.assign(3000,0.);
    for(int y=1;y<=h->GetNbinsY();++y)for(int x=1;x<=3000;++x) {
      double n=h->GetBinContent(x,y);require(std::isfinite(n) && n>=0.,"Invalid histogram count");
      (y%2?s.fit:s.check)[x-1]+=n;
    }
  }
  require(std::all_of(found.begin(),found.end(),[](bool b){return b;}),"Incomplete channels");
  return out;
}
struct Cdf {
  std::vector<double> sum;
  explicit Cdf(const std::vector<double>& bins) {
    // A fixed 50-keV kernel removes numerical discontinuities; conserve counts.
    std::vector<double> smooth(bins.size(),0.);
    for(size_t i=0;i<bins.size();++i) {
      smooth[i]+=.5*bins[i];
      smooth[i?i-1:i]+=.25*bins[i];smooth[i+1<bins.size()?i+1:i]+=.25*bins[i];
    }
    sum.assign(bins.size()+1,0.);std::partial_sum(smooth.begin(),smooth.end(),sum.begin()+1);
  }
  double operator()(double energy) const {
    if(energy<=0.)return 0.;
    double x=energy/bin_width;
    if(x>=sum.size()-1)return sum.back();
    size_t i=size_t(x);return sum[i]+(x-i)*(sum[i+1]-sum[i]);
  }
};
struct Metric {
  const Cdf& target;
  std::vector<double> root_p;
  double lo,hi,ref_count=0.;
  Metric(const Cdf& reference,const Cdf& current,double low=fit_low,double high=fit_high):target(current),lo(low),hi(high) {
    int n=int(std::llround((hi-lo)/fit_bin));root_p.resize(n);
    for(int i=0;i<n;++i){root_p[i]=std::max(0.,reference(lo+(i+1)*fit_bin)-reference(lo+i*fit_bin));ref_count+=root_p[i];}
    if(ref_count>0.)for(auto& v:root_p)v=std::sqrt(v/ref_count);
  }
  double count(double k) const {return std::max(0.,target(hi/k)-target(lo/k));}
  double operator()(double k) const {
    const double n=count(k);if(ref_count<1. || n<minimum_count/4.)return 10.;
    double affinity=0.,last=target(lo/k);
    for(size_t i=0;i<root_p.size();++i) {
      double next=target((lo+(i+1)*fit_bin)/k);
      affinity+=root_p[i]*std::sqrt(std::max(0.,next-last)/n);last=next;
    }
    return std::max(0.,1.-affinity);
  }
};
struct Estimate {double k=1.,distance=10.,support_low=.5,support_high=5.;bool boundary=false;std::vector<double> grid,values;};
Estimate estimate(const Metric& metric) {
  Estimate e;int best=0;
  for(int i=0;i<points;++i) {
    double k=std::exp(std::log(.5)+i*std::log(10.)/(points-1));
    e.grid.push_back(k);e.values.push_back(metric(k));if(e.values[i]<e.values[best])best=i;
  }
  e.boundary=best==0 || best==points-1;
  double lo=e.grid[std::max(0,best-1)],hi=e.grid[std::min(points-1,best+1)];
  const double ratio=(std::sqrt(5.)-1.)/2.;
  double c=hi-ratio*(hi-lo),d=lo+ratio*(hi-lo),fc=metric(c),fd=metric(d);
  for(int i=0;i<36;++i) {
    if(fc<fd){hi=d;d=c;fd=fc;c=hi-ratio*(hi-lo);fc=metric(c);}
    else{lo=c;c=d;fc=fd;d=lo+ratio*(hi-lo);fd=metric(d);}
  }
  e.k=(lo+hi)/2.;e.distance=metric(e.k);
  std::vector<double> support;
  for(int i=0;i<points;++i)if(e.values[i]<=e.distance+.005)support.push_back(e.grid[i]);
  if(!support.empty()){e.support_low=support.front();e.support_high=support.back();}
  return e;
}
}
int main(int argc,char** argv) {
  try {
    fs::path input,output,run_list;
    for(int i=1;i<argc;++i) {
      std::string a=argv[i];
      if(a=="--help"){std::cout<<"calibrate --input-dir DIR --output-dir NEW_DIR --run-list FILE\n";return 0;}
      auto value=[&](){require(++i<argc,"Missing option value");return std::string(argv[i]);};
      if(a=="--input-dir")input=value();else if(a=="--output-dir")output=value();else if(a=="--run-list")run_list=value();
      else throw std::runtime_error("Unknown argument "+a);
    }
    require(!input.empty() && !output.empty() && !run_list.empty(),"Use --help");
    require(!fs::exists(output),"Output directory already exists");fs::create_directories(output);
    std::ifstream list(run_list);require(bool(list),"Cannot open run list");
    std::vector<int> runs;int number;while(list>>number)runs.push_back(number);
    require(list.eof() && !runs.empty(),"Invalid run list");std::sort(runs.begin(),runs.end());
    require(std::adjacent_find(runs.begin(),runs.end())==runs.end(),"Duplicate run");
    const std::map<std::string,int> reference_runs={{"high1",472},{"low",616},{"high2",680}};
    std::map<std::string,Input> references;
    for(const auto& kv:reference_runs)references.emplace(kv.first,load(input/("spectra_"+rn(kv.second)+".root")));
    json rows=json::array(),sources=json::array();
    std::unique_ptr<TFile> diagnostics(TFile::Open((output/"profiles.root").c_str(),"CREATE","Descriptive distance profiles",404));
    require(diagnostics && !diagnostics->IsZombie(),"Cannot create diagnostics");
    for(int run:runs) {
      auto data=load(input/("spectra_"+rn(run)+".root"));std::string period=data.meta.at("period");
      auto& reference=references.at(period);int ref_run=reference_runs.at(period);
      sources.push_back({{"run",rn(run)},{"path",(input/("spectra_"+rn(run)+".root")).string()},{"sha256",data.sha},
        {"sample_entries",data.meta.at("sample_entries")},{"input_entries",data.meta.at("input_entries")},{"sample_fraction",data.meta.at("sample_fraction")}});
      for(int i=0;i<channels;++i) {
        const auto& current=data.spectra[i];const auto& ref=reference.spectra[i];
        require(current.name==ref.name && ref.old_k==1.,"Inconsistent channel/reference scale");
        Cdf rt(ref.fit),rv(ref.check),tt(current.fit),tv(current.check);
        Metric train(rt,tt),validation(rv,tv),ref_self(rt,rv),alternate(rt,tt,30000.,80000.);
        auto fit=estimate(train),check=estimate(validation),window=estimate(alternate);
        double old_fit=train(current.old_k),old_val=validation(current.old_k),new_val=validation(fit.k);
        double split_pct=100.*std::abs(std::log(check.k/fit.k));
        double window_pct=100.*std::abs(std::log(window.k/fit.k));
        double support_pct=100.*std::log(fit.support_high/fit.support_low);
        std::vector<std::string> flags;
        if(train.ref_count<minimum_count || validation.ref_count<minimum_count)flags.push_back("insufficient_reference");
        if(train.count(fit.k)<minimum_count || validation.count(fit.k)<minimum_count)flags.push_back("insufficient_target");
        if(fit.boundary || check.boundary)flags.push_back("search_boundary");
        if(support_pct>15.)flags.push_back("weak_scale_identification");
        if(split_pct>3.)flags.push_back("fit_check_scale_disagreement");
        if(window_pct>3.)flags.push_back("energy_window_sensitivity");
        if(new_val>std::max(.08,3.*ref_self(1.)+.015))flags.push_back("shape_mismatch");
        bool improve=new_val<old_val*.9 && old_val-new_val>.002;
        std::string decision=flags.empty()?(improve?"improved_candidate":"keep_previous_no_clear_gain"):"keep_previous_flagged";
        double proposed=(flags.empty() && improve)?fit.k:current.old_k;
        if(run==ref_run){fit.k=1.;proposed=1.;new_val=validation(1.);decision="fixed_reference";flags.clear();}
        rows.push_back({{"run",rn(run)},{"period",period},{"reference_run",rn(ref_run)},{"channel",current.name},{"index",i},
          {"old_k",current.old_k},{"candidate_k",fit.k},{"proposed_k",proposed},{"decision",decision},{"flags",flags},
          {"fit_distance_old",old_fit},{"fit_distance_candidate",train(fit.k)},{"check_distance_old",old_val},{"check_distance_candidate",new_val},
          {"reference_internal_distance",ref_self(1.)},{"reference_fit_count",train.ref_count},{"reference_check_count",validation.ref_count},
          {"target_fit_count",train.count(fit.k)},{"target_check_count",validation.count(fit.k)},
          {"independent_check_k",check.k},{"alternate_window_k",window.k},{"split_difference_pct",split_pct},{"window_difference_pct",window_pct},
          {"profile_support_k",{fit.support_low,fit.support_high}},{"profile_support_log_width_pct",support_pct},
          {"uncertainty_status","No calibrated confidence interval; split/window/profile spread are descriptive diagnostics"}});
        diagnostics->cd();TGraph profile(points,fit.grid.data(),fit.values.data());
        profile.SetName(("profile_"+rn(run)+"_"+current.name).c_str());profile.SetTitle("Fit-only Hellinger distance;k;H squared");profile.Write();
      }
      std::cout<<"complete run="<<rn(run)<<" rows="<<rows.size()<<std::endl;
    }
    json settings={{"method","Normalized spectrum registration: minimize squared Hellinger distance"},
      {"formula","H2 = 1 - sum(sqrt(p_i*q_i(k))); p and q separately normalized inside the comparison window"},
      {"fit_window_keV",{fit_low,fit_high}},{"check_window_keV",{fit_low,fit_high}},{"alternate_fit_window_keV",{30000,80000}},
      {"comparison_bin_keV",fit_bin},{"source_bin_keV",bin_width},{"smoothing","[0.25,0.5,0.25] at 50 keV, counts conserved"},
      {"k_search_range",{.5,5.}},{"grid_points",points},{"minimum_window_counts",minimum_count},
      {"maximum_split_log_difference_pct",3.},{"maximum_window_log_difference_pct",3.},{"maximum_profile_support_log_width_pct",15.},
      {"support_distance_increment",.005},{"minimum_relative_validation_gain",.10},{"minimum_absolute_validation_gain",.002},
      {"shape_mismatch_limit","max(0.08,3*reference_internal_H2+0.015)"},
      {"reference_runs",reference_runs},{"adoption","Candidates only; no production parameter changes"},
      {"limitations","Conditional shape distance, not likelihood or chi-square; changing k changes window membership. Shape changes can invalidate a pure scale model."}};
    json result={{"status","complete"},{"settings",settings},{"sources",sources},{"rows",rows},
      {"program_sha256",hash(read_text(fs::canonical("/proc/self/exe")))}};
    TNamed provenance("metadata",settings.dump().c_str());provenance.Write();diagnostics->Close();
    std::ofstream report(output/"results.json");report<<result.dump(2)<<'\n';report.close();require(bool(report),"Cannot save results");
    std::ofstream csv(output/"parameters.csv");csv<<"run,period,reference_run,channel,old_k,candidate_k,proposed_k,decision,flags\n";csv<<std::setprecision(17);
    fs::create_directories(output/"candidate_parameters");
    for(int i=0;i<channels;++i) {
      std::string name=i<24?"ring"+rn(i+1).substr(2):"sector"+rn(i-23).substr(2);
      std::ofstream table(output/"candidate_parameters"/(name+".txt"));table<<std::setprecision(17);
      for(const auto& r:rows)if(int(r.at("index"))==i)table<<std::stoi(r.at("run").get<std::string>())<<"  "<<r.at("proposed_k").get<double>()<<'\n';
      table.close();require(bool(table),"Cannot export candidate table");
    }
    for(const auto& r:rows) {
      csv<<r.at("run").get<std::string>()<<','<<r.at("period").get<std::string>()<<','<<r.at("reference_run").get<std::string>()<<','
        <<r.at("channel").get<std::string>()<<','<<r.at("old_k").get<double>()<<','<<r.at("candidate_k").get<double>()<<','
        <<r.at("proposed_k").get<double>()<<','<<r.at("decision").get<std::string>()<<',';
      for(size_t j=0;j<r.at("flags").size();++j)csv<<(j?"|":"")<<r.at("flags")[j].get<std::string>();
      csv<<'\n';
    }
    csv.close();require(bool(csv),"Cannot save CSV");return 0;
  }catch(const std::exception& e){std::cerr<<"s3 attenuation registration: "<<e.what()<<'\n';return 1;}
}
