// Fit and check S3 normalization on independent source-entry blocks.
// The extractor retains pre-normalization energies and raw multiplicities.
#include "sm150_event_calibration.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TObjString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <array>
#include <chrono>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>

using namespace sm150::diagnostic;
namespace {
struct Options {
  fs::path input,outdir;
  double low=15000.,high=100000.,time_gate=80.,min_span=5000.;
  size_t min_points=100;
};
struct Point {double x,y;};
struct Fit {
  bool ok=false;std::string reason="not_fitted";size_t count=0,inliers=0;
  int iterations=0;double a=0.,b=1.,scale=0.,span=0.;
  std::array<double,3> covariance{}; // Var(a), Cov(a,b), Var(b), conditional on weights.
};
double quantile(std::vector<double> v,double p) {
  if(v.empty())return std::numeric_limits<double>::quiet_NaN();
  const double index=p*(v.size()-1);const size_t lo=size_t(index),hi=std::min(lo+1,v.size()-1);
  std::nth_element(v.begin(),v.begin()+lo,v.end());const double a=v[lo];
  if(hi==lo)return a;
  std::nth_element(v.begin(),v.begin()+hi,v.end());
  return a+(index-lo)*(v[hi]-a);
}
Fit robust_line(const std::vector<Point>& points,const Options& o) {
  Fit f;f.count=points.size();
  if(f.count<o.min_points){f.reason="insufficient_training_points";return f;}
  std::vector<double> x;for(const auto& p:points)x.push_back(p.x);
  f.span=quantile(x,.9)-quantile(x,.1);
  if(f.span<o.min_span){f.reason="insufficient_energy_span";return f;}
  // Binned medians give a deterministic start, independent of old coefficients.
  std::map<int,std::pair<std::vector<double>,std::vector<double>>> bins;
  for(const auto& p:points){auto& b=bins[int(std::floor(p.x/2500.))];b.first.push_back(p.x);b.second.push_back(p.y);}
  std::vector<Point> medians;
  for(const auto& b:bins)if(b.second.first.size()>=5)medians.push_back({quantile(b.second.first,.5),quantile(b.second.second,.5)});
  if(medians.size()<3){f.reason="insufficient_populated_energy_bins";return f;}
  std::vector<double> slopes;
  for(size_t i=0;i<medians.size();++i)for(size_t j=i+1;j<medians.size();++j)
    if(medians[j].x-medians[i].x>1000.)slopes.push_back((medians[j].y-medians[i].y)/(medians[j].x-medians[i].x));
  if(slopes.empty()){f.reason="singular_start";return f;}
  f.b=quantile(slopes,.5);std::vector<double> intercepts;
  for(const auto& p:points)intercepts.push_back(p.y-f.b*p.x);
  f.a=quantile(intercepts,.5);
  bool converged=false;
  double sumw=0.,meanx=0.,sxx=0.,weighted_rss=0.;
  for(int iteration=0;iteration<50;++iteration) {
    std::vector<double> residual;for(const auto& p:points)residual.push_back(p.y-f.a-f.b*p.x);
    const double center=quantile(residual,.5);
    std::vector<double> absolute;for(double r:residual)absolute.push_back(std::abs(r-center));
    f.scale=std::max(5.,1.4826*quantile(absolute,.5));
    sumw=0.;double sx=0.,sy=0.;std::vector<double> weights;weights.reserve(points.size());
    for(size_t i=0;i<points.size();++i) {
      double w=std::abs(residual[i])<=1.5*f.scale?1.:1.5*f.scale/std::abs(residual[i]);
      weights.push_back(w);sumw+=w;sx+=w*points[i].x;sy+=w*points[i].y;
    }
    if(sumw<=2){f.reason="insufficient_effective_weight";return f;}
    meanx=sx/sumw;const double meany=sy/sumw;sxx=0.;double sxy=0.;
    for(size_t i=0;i<points.size();++i) {
      const double dx=points[i].x-meanx;sxx+=weights[i]*dx*dx;sxy+=weights[i]*dx*(points[i].y-meany);
    }
    if(sxx<=0 || !std::isfinite(sxx)){f.reason="singular_design";return f;}
    const double b=sxy/sxx,a=meany-b*meanx;
    const double shift=std::max(std::abs(a-f.a+(b-f.b)*o.low),std::abs(a-f.a+(b-f.b)*o.high));
    f.a=a;f.b=b;f.iterations=iteration+1;
    weighted_rss=0.;f.inliers=0;
    for(size_t i=0;i<points.size();++i) {
      const double r=points[i].y-a-b*points[i].x;weighted_rss+=weights[i]*r*r;
      if(std::abs(r)<=2.5*f.scale)++f.inliers;
    }
    if(shift<1e-3){converged=true;break;}
  }
  if(!converged){f.reason="iteration_limit";return f;}
  if(!std::isfinite(f.a) || !std::isfinite(f.b) || f.b<=0){f.reason="nonpositive_or_nonfinite_slope";return f;}
  const double variance=weighted_rss/(sumw-2.);
  f.covariance={variance*(1./sumw+meanx*meanx/sxx),-variance*meanx/sxx,variance/sxx};
  f.ok=true;f.reason="fit";return f;
}
struct Model {
  int face=0,id=0,key=0;double old_a=0,old_b=1,a=0,b=1;
  std::array<double,3> covariance{};bool measured=false,fixed=false;
  std::string status="inherited_not_fitted";Fit local;
};
Model compose(Model m,const Fit& f,const Model& reference) {
  m.local=f;
  if(!f.ok){m.status="inherited_"+f.reason;return m;}
  m.a=reference.a+reference.b*f.a;m.b=reference.b*f.b;
  const auto& c=reference.covariance;const auto& q=f.covariance;
  m.covariance={c[0]+f.a*f.a*c[2]+2*f.a*c[1]+reference.b*reference.b*q[0],
                f.b*c[1]+f.a*f.b*c[2]+reference.b*reference.b*q[1],
                f.b*f.b*c[2]+reference.b*reference.b*q[2]};
  m.measured=true;m.status="candidate_fit";return m;
}
json fit_json(const Fit& f) {
  return {{"status",f.reason},{"n_train",f.count},{"n_core",f.inliers},{"iterations",f.iterations},
          {"x_10_90_span_keV",f.span},{"residual_mad_sigma_keV",f.scale},
          {"local_p0_keV",f.ok?json(f.a):json(nullptr)},{"local_p1",f.ok?json(f.b):json(nullptr)}};
}
json metrics(const std::vector<double>& v) {
  if(v.empty())return {{"n",0},{"median_keV",nullptr},{"half_central68_keV",nullptr},{"fraction_outside_2500keV",nullptr}};
  size_t outside=0;for(double x:v)if(std::abs(x)>2500)++outside;
  return {{"n",v.size()},{"median_keV",quantile(v,.5)},
          {"half_central68_keV",.5*(quantile(v,.84)-quantile(v,.16))},
          {"fraction_outside_2500keV",double(outside)/v.size()}};
}
std::string number(int n){std::ostringstream s;s<<std::setw(2)<<std::setfill('0')<<n;return s.str();}
std::string file_sha(const fs::path& p) {
  std::ifstream in(p,std::ios::binary);require(bool(in),"Cannot hash "+p.string());
  Digest d;std::vector<char> data(1024*1024);
  while(in){in.read(data.data(),data.size());const auto n=in.gcount();if(n)d.add(data.data(),n);}
  require(in.eof(),"Input hashing failed");return d.finish();
}
struct ChannelHist {
  std::unique_ptr<TH2D> fit;
  std::array<std::unique_ptr<TH2D>,3> residual;
  std::array<std::unique_ptr<TH1D>,3> one;
  std::array<std::vector<double>,3> values,reference_values;
  std::vector<std::array<std::vector<double>,2>> blocks;
  ChannelHist(const std::string& name,int nblocks,const Options& o) :blocks(nblocks) {
    fit=std::make_unique<TH2D>((name+"_fit").c_str(),"Training reference pairs;local pre-fb energy (keV);reference pre-fb energy (keV)",170,o.low,o.high,200,0,o.high);
    fit->SetDirectory(nullptr);
    const std::string state[]={"before","previous","candidate"};
    for(int k=0;k<3;++k) {
      residual[k]=std::make_unique<TH2D>((name+"_"+state[k]+"_2d").c_str(),"Independent validation;local pre-fb energy (keV);sector-ring energy (keV)",85,o.low,o.high,2000,-20000,20000);
      one[k]=std::make_unique<TH1D>((name+"_"+state[k]+"_1d").c_str(),"Held-out validation;sector-ring energy (keV);counts / 20 keV",10000,-100000,100000);
      residual[k]->SetDirectory(nullptr);one[k]->SetDirectory(nullptr);
    }
  }
};
}

int main(int argc,char** argv) {
  try {
    Options o;
    for(int i=1;i<argc;++i) {
      const std::string a=argv[i];
      if(a=="--help") {
        std::cout<<"./calibrate --input s3.root --outdir NEW_DIR [--min-points 100] [--min-span-kev 5000]\n"
          "[--energy-low-kev 15000 --energy-high-kev 100000 --time-gate-ns 80]\n"
          "Even blocks fit; odd blocks validate. Candidate outputs never replace production parameters.\n";return 0;
      }
      require(i+1<argc,"Missing option value");std::string v=argv[++i];
      if(a=="--input")o.input=v;else if(a=="--outdir")o.outdir=v;
      else if(a=="--min-points")o.min_points=std::stoull(v);
      else if(a=="--min-span-kev")o.min_span=std::stod(v);
      else if(a=="--energy-low-kev")o.low=std::stod(v);
      else if(a=="--energy-high-kev")o.high=std::stod(v);
      else if(a=="--time-gate-ns")o.time_gate=std::stod(v);
      else throw std::runtime_error("Unknown argument: "+a);
    }
    require(!o.input.empty() && !o.outdir.empty() && o.min_points>=3 && std::isfinite(o.min_span) && o.min_span>0 && std::isfinite(o.low) && std::isfinite(o.high) && std::isfinite(o.time_gate) && o.low>=0 && o.high>o.low && o.time_gate>0,"Invalid fitting settings");
    o.input=fs::canonical(o.input);o.outdir=fs::absolute(o.outdir);
    require(!fs::exists(o.outdir),"Output directory already exists; use a new candidate directory");
    const Stat before(o.input);const std::string input_hash=file_sha(o.input);
    std::unique_ptr<TFile> input(TFile::Open(o.input.c_str(),"READ"));require(input && !input->IsZombie(),"Cannot read S3 calibration ROOT");
    auto* object=dynamic_cast<TObjString*>(input->Get("metadata"));require(object,"Missing extraction provenance");
    const json meta=json::parse(object->GetString().Data());
    require(meta.at("schema")=="sm150_s3_calibration_v1" && meta.at("status")=="complete" && meta.at("energy_agreement_gate")==false && meta.at("addback")==false,"Unsupported calibration source");
    require(o.time_gate<=meta.at("cluster_gap_ns").get<double>(),"Fit time gate exceeds extraction coverage");
    auto* tree=dynamic_cast<TTree*>(input->Get("s3_pairs"));require(tree,"Missing S3 pairs");
    std::array<Model,56> models;
    std::array<bool,56> seen{};
    for(auto it=meta.at("channels").begin();it!=meta.at("channels").end();++it) {
      const auto& c=it.value();const int kind=c.at("kind");if(kind!=3 && kind!=4)continue;
      const int face=kind-3,id=c.at(face?"sector":"ring"),index=face?24+id-1:id-1;
      require(index>=0 && index<56 && !seen[index],"Invalid channel mapping");seen[index]=true;
      auto& m=models[index];m.face=face;m.id=id;m.key=std::stoi(it.key());
      m.a=m.old_a=c.at("fb_p0_kev");m.b=m.old_b=c.at("fb_p1");
    }
    require(std::all_of(seen.begin(),seen.end(),[](bool v){return v;}),"Missing S3 calibration channels");
    auto& anchor=models[24];anchor.fixed=true;anchor.status="fixed_existing_sector1_reference";
    TTreeReader reader(tree);
    TTreeReaderValue<Int_t> ring(reader,"ring_id"),sector(reader,"sector_id"),block(reader,"block"),flags(reader,"quality_flags");
    TTreeReaderValue<Bool_t> boundary(reader,"boundary");
    TTreeReaderValue<Double_t> r(reader,"ring_energy_pre_fb_keV"),s(reader,"sector_energy_pre_fb_keV"),
      ro(reader,"ring_energy_previous_keV"),so(reader,"sector_energy_previous_keV"),dt(reader,"dt_sector_minus_ring_ns");
    uint64_t read=0,selected=0;std::array<std::vector<Point>,56> train;
    const int blocks=meta.at("blocks");
    while(reader.Next()) {
      ++read;require(*ring>=1 && *ring<=24 && *sector>=1 && *sector<=32 && *block>=0 && *block<blocks,"Invalid pair identity");
      require(std::isfinite(*r) && std::isfinite(*s) && std::isfinite(*ro) && std::isfinite(*so) && std::isfinite(*dt),"Nonfinite pair");
      if(*flags || *boundary || std::abs(*dt)>o.time_gate || *r<o.low || *r>=o.high || *s<o.low || *s>=o.high)continue;
      ++selected;
      if(*block%2==0) {
        if(*sector==1)train[*ring-1].push_back({*r,*s});
        if(*ring==11 && *sector!=1)train[24+*sector-1].push_back({*s,*r});
      }
    }
    require(read==uint64_t(tree->GetEntries()) && read==meta.at("single_pairs").get<uint64_t>(),"Input pair count mismatch");
    for(int i=0;i<24;++i)models[i]=compose(models[i],robust_line(train[i],o),anchor);
    const Model ring_reference=models[10];
    for(int i=25;i<56;++i) {
      if(!ring_reference.measured)models[i].status="inherited_ring11_reference_not_measured";
      else models[i]=compose(models[i],robust_line(train[i],o),ring_reference);
    }
    // Freeze every coefficient before examining the held-out block residuals.
    fs::create_directories(o.outdir);
    std::unique_ptr<TFile> diagnostics(TFile::Open((o.outdir/"diagnostics.root").c_str(),"CREATE"));
    require(diagnostics && !diagnostics->IsZombie(),"Cannot create diagnostics");
    std::array<std::unique_ptr<ChannelHist>,56> hist;
    for(int i=0;i<56;++i)hist[i]=std::make_unique<ChannelHist>(std::string(i<24?"ring":"sector")+number(models[i].id),blocks,o);
    const std::string states[]={"before","previous","candidate"};
    std::array<std::unique_ptr<TH2D>,3> global;
    std::array<std::unique_ptr<TH1D>,3> global_residual;
    std::array<std::vector<double>,3> global_values;
    for(int k=0;k<3;++k) {
      global[k]=std::make_unique<TH2D>(("global_"+states[k]).c_str(),"Held-out single-hit pairs;ring energy (keV);sector energy (keV)",200,0,100000,200,0,100000);
      global_residual[k]=std::make_unique<TH1D>(("global_"+states[k]+"_residual").c_str(),"Held-out pairs;sector-ring energy (keV);counts / 20 keV",10000,-100000,100000);
      global[k]->SetDirectory(nullptr);global_residual[k]->SetDirectory(nullptr);
    }
    for(int i=0;i<56;++i)for(const auto& p:train[i])hist[i]->fit->Fill(p.x,p.y);
    uint64_t n_validation=0,n_train=0;
    // A second sequential pass avoids holding every selected pair in memory.
    // Exact residual quantiles are retained; no random downsampling is used.
    reader.Restart();
    uint64_t read_validation=0;
    while(reader.Next()) {
      ++read_validation;
      if(*flags || *boundary || std::abs(*dt)>o.time_gate || *r<o.low || *r>=o.high || *s<o.low || *s>=o.high)continue;
      if(*block%2==0){++n_train;continue;}++n_validation;
      const int ri=*ring-1,si=24+*sector-1;
      const auto& rm=models[ri];const auto& sm=models[si];
      const double er[]={*r,*ro,rm.a+rm.b*(*r)},es[]={*s,*so,sm.a+sm.b*(*s)};
      for(int k=0;k<3;++k) {
        const double difference=es[k]-er[k];global[k]->Fill(er[k],es[k]);global_residual[k]->Fill(difference);global_values[k].push_back(difference);
        for(int index:{ri,si}) {
          hist[index]->residual[k]->Fill(index<24?*r:*s,difference);
          hist[index]->one[k]->Fill(difference);hist[index]->values[k].push_back(difference);
          if((index<24 && *sector==1) || (index>=24 && *ring==11))hist[index]->reference_values[k].push_back(difference);
          if(k)hist[index]->blocks[*block][k-1].push_back(difference);
        }
      }
    }
    require(read_validation==read && n_train+n_validation==selected,"Validation pass count mismatch");
    json rows=json::array(),per_block=json::array();int fitted=0,inherited=0;
    for(int i=0;i<56;++i) {
      const auto& m=models[i];if(m.measured)++fitted;else if(!m.fixed)++inherited;
      json row={{"face",m.face?"sector":"ring"},{"id",m.id},{"channel",m.key},{"status",m.status},
        {"old_p0_keV",m.old_a},{"old_p1",m.old_b},{"candidate_p0_keV",m.a},{"candidate_p1",m.b},
        {"fit",fit_json(m.local)},{"reference",m.fixed?"fixed existing sector1 scale":m.face?"corrected ring11":"fixed existing sector1"},
        {"uncertainty_kind",m.fixed?"fixed reference; absolute scale uncertainty not estimated":m.measured?"conditional weighted regression; shared ring11 reference propagated":"unmeasured/inherited; old uncertainty unavailable"},
        {"conditional_covariance",m.measured?json(m.covariance):json(nullptr)},
        {"p0_error_keV",m.measured?json(std::sqrt(std::max(0.,m.covariance[0]))):json(nullptr)},
        {"p1_error",m.measured?json(std::sqrt(std::max(0.,m.covariance[2]))):json(nullptr)}};
      for(int k=0;k<3;++k) {
        row["validation"][states[k]]=metrics(hist[i]->values[k]);
        row["validation_reference"][states[k]]=metrics(hist[i]->reference_values[k]);
      }
      const auto& old=hist[i]->values[1];const auto& now=hist[i]->values[2];
      row["validation_sufficient"]=old.size()>=50;
      if(old.size()>=50) {
        const double om=std::abs(quantile(old,.5)),nm=std::abs(quantile(now,.5));
        const double ow=.5*(quantile(old,.84)-quantile(old,.16)),nw=.5*(quantile(now,.84)-quantile(now,.16));
        row["validation_assessment"]=(nm<=om+100. && nw<=1.1*ow+50.)?"no_large_degradation_in_this_split":"review_degradation";
      }else row["validation_assessment"]="insufficient_validation_events";
      rows.push_back(row);
      for(int b=1;b<blocks;b+=2)per_block.push_back({{"face",m.face?"sector":"ring"},{"id",m.id},{"block",b},
        {"previous",metrics(hist[i]->blocks[b][0])},{"candidate",metrics(hist[i]->blocks[b][1])}});
      diagnostics->cd();hist[i]->fit->Write();for(int k=0;k<3;++k){hist[i]->residual[k]->Write();hist[i]->one[k]->Write();}
    }
    for(int k=0;k<3;++k){global[k]->Write();global_residual[k]->Write();}
    json summary={{"status","complete"},{"run",meta.at("run")},{"period",meta.at("period")},
      {"source",o.input.string()},{"source_sha256",input_hash},{"raw_source",meta.at("input")},
      {"raw_entries_read",meta.at("processed_entries")},{"full_raw_input_read",meta.at("full_input_read")},
      {"timing_profile",meta.at("timing_profile")},{"timing_source_sha256",meta.at("timing_source_sha256")},
      {"input_pair_rows",read},{"selected_pairs",selected},{"fit_block_pairs",n_train},{"validation_pairs",n_validation},
      {"fit_energy_pre_fb_keV",{o.low,o.high}},{"time_gate_abs_ns",o.time_gate},
      {"min_training_points",o.min_points},{"min_10_90_span_keV",o.min_span},
      {"reference","Existing sector1 energy scale fixed; ring11 carries the second-stage reference"},
      {"fit_method","Binned-median Theil-Sen start; Huber IRLS, k=1.5, MAD scale floor 5 keV, max50 iterations"},
      {"covariance_scope","Approximate conditional weighted LS; ring11 statistical reference propagated to sectors. No absolute calibration, robust selection or temporal/model systematic uncertainty included."},
      {"selection","Unmerged raw multiplicity1+1, clean flags, nonboundary, fixed pre-fb energy/time cuts. No energy-agreement gate. Even blocks fit; odd blocks only validate."},
      {"candidate_fitted_channels",fitted},{"inherited_channels",inherited},{"fixed_reference_channels",1},
      {"parameters",rows},{"block_validation",per_block},{"production_parameters_modified",false}};
    for(int k=0;k<3;++k)summary["validation_global"][states[k]]=metrics(global_values[k]);
    const auto encoded=summary.dump();TObjString(encoded.c_str()).Write("summary");diagnostics->Close();
    for(int face=0;face<2;++face) {
      const auto name=std::string("correction_")+(face?"sector_":"ring_")+meta.at("run").get<std::string>()+".txt";
      std::ofstream table(o.outdir/name);table<<"# candidate only: id p0_keV p1; consult results.json for fit/inherited status\n"<<std::setprecision(17);
      for(const auto& m:models)if(m.face==face)table<<m.id<<' '<<m.a<<' '<<m.b<<'\n';
      table.close();require(bool(table),"Cannot export candidate table");
    }
    require(before.same(Stat(o.input)),"Calibration source changed during fitting");
    std::ofstream report(o.outdir/"results.json");report<<summary.dump(2)<<'\n';report.close();require(bool(report),"Cannot save fit results");
    std::cout<<"complete run="<<meta.at("run")<<" selected="<<selected<<" fitted="<<fitted<<" inherited="<<inherited<<std::endl;
    return 0;
  }catch(const std::exception& e){std::cerr<<"s3_fb_cor: "<<e.what()<<'\n';return 1;}
}
