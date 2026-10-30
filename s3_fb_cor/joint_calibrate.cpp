// Robust pixel correlations followed by a connected global strip calibration.
// The common energy scale and zero are fixed explicitly, not inferred from Efront=Eback.
#include "sm150_event_calibration.h"
#include <TFile.h>
#include <TObjString.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <queue>
#include <unistd.h>
using namespace sm150::diagnostic;
namespace {
struct Point { double x=0,y=0; int block=0; };
struct Edge {
  int r=0,s=0,n=0,random=0,iterations=0;
  bool ok=false;std::string status="not_fit";
  double a=0,b=1,scale=0,span=0,xlo=0,xhi=0,neff=0;
};
double quantile(std::vector<double> v,double p) {
  if(v.empty())return std::numeric_limits<double>::quiet_NaN();
  double pos=p*(v.size()-1);size_t lo=size_t(pos),hi=std::min(lo+1,v.size()-1);
  std::nth_element(v.begin(),v.begin()+lo,v.end());double a=v[lo];
  std::nth_element(v.begin(),v.begin()+hi,v.end());return a+(v[hi]-a)*(pos-lo);
}
Edge fit_edge(const std::vector<Point>& points,int r,int s,int random,int subset,int minimum) {
  Edge e;e.r=r;e.s=s;e.random=random;std::vector<Point> p;std::vector<double> xs,ys;
  for(const auto& v:points)if(subset<0 || v.block%4==subset){p.push_back(v);xs.push_back(v.x);ys.push_back(v.y);}
  e.n=int(p.size());if(e.n<minimum){e.status="too_few_pairs";return e;}
  e.xlo=quantile(xs,.10);e.xhi=quantile(xs,.90);e.span=e.xhi-e.xlo;
  if(e.span<5.){e.status="energy_span_below_5MeV";return e;}
  if(subset<0 && e.n-.5*random<.5*e.n){e.status="large_time_sideband_fraction";return e;}
  // Median energy slices initialize the main correlation before robust TLS.
  std::map<int,std::vector<Point>> bins;
  for(const auto& v:p)bins[int(v.x/2.5)].push_back(v);
  std::vector<Point> medians;
  for(const auto& kv:bins)if(kv.second.size()>=8) {
    std::vector<double> a,b;for(const auto& v:kv.second){a.push_back(v.x);b.push_back(v.y);}
    medians.push_back({quantile(a,.5),quantile(b,.5),0});
  }
  std::vector<double> slopes;
  for(size_t i=0;i<medians.size();++i)for(size_t j=i+1;j<medians.size();++j)
    if(medians[j].x-medians[i].x>=5.)slopes.push_back((medians[j].y-medians[i].y)/(medians[j].x-medians[i].x));
  e.b=slopes.empty()?1.:quantile(slopes,.5);
  if(!(e.b>.1 && e.b<10.))e.b=1.;
  std::vector<double> offsets;for(const auto& v:p)offsets.push_back(v.y-e.b*v.x);e.a=quantile(offsets,.5);
  double sumw=0,sumw2=0;
  for(int iteration=0;iteration<60;++iteration) {
    std::vector<double> residual;const double norm=std::hypot(1.,e.b);
    for(const auto& v:p)residual.push_back((v.y-e.a-e.b*v.x)/norm);
    const double center=quantile(residual,.5);std::vector<double> deviation;
    for(double d:residual)deviation.push_back(std::abs(d-center));
    const double scale=std::max(.015,1.4826*quantile(deviation,.5));
    sumw=sumw2=0;double sx=0,sy=0,sxx=0,syy=0,sxy=0;
    for(size_t i=0;i<p.size();++i) {
      double z=residual[i]/(2.385*scale),w=1./(1.+z*z);const auto& v=p[i];
      sumw+=w;sumw2+=w*w;sx+=w*v.x;sy+=w*v.y;sxx+=w*v.x*v.x;syy+=w*v.y*v.y;sxy+=w*v.x*v.y;
    }
    double vx=sxx/sumw-std::pow(sx/sumw,2),vy=syy/sumw-std::pow(sy/sumw,2),cov=sxy/sumw-sx*sy/(sumw*sumw);
    if(cov<=0 || vx<=0 || vy<=0){e.status="nonpositive_correlation";return e;}
    double b=(vy-vx+std::hypot(vy-vx,2*cov))/(2*cov),a=(sy-b*sx)/sumw;
    double change=std::abs(a-e.a)+50*std::abs(b-e.b);e.a=a;e.b=b;e.scale=scale;e.iterations=iteration+1;
    if(change<1e-6)break;
  }
  e.neff=sumw*sumw/sumw2;
  if(!(e.b>.2 && e.b<5.)){e.status="slope_out_of_range";return e;}
  if(e.scale>.3*e.span){e.status="unresolved_wide_cloud";return e;}
  e.ok=true;e.status="fit";return e;
}
json edge_json(const Edge& e) {
  return {{"r",e.r},{"s",e.s},{"n",e.n},{"random",e.random},{"ok",e.ok},{"status",e.status},{"a_MeV",e.a},{"b",e.b},
    {"orthogonal_scale_MeV",e.scale},{"span_MeV",e.span},{"xlo_MeV",e.xlo},{"xhi_MeV",e.xhi},{"neff",e.neff},{"iterations",e.iterations}};
}
struct Model {
  std::array<double,56> a{},b{};std::array<std::string,56> status{};
  double condition=0;int iterations=0,rank=0;bool ok=false;
  Model(){b.fill(1.);status.fill("unchanged_disconnected");}
};
struct TemporalPrior {int index=0;double a=0,b=1,error=.5;};
Model global_fit(const std::vector<Edge>& edges,int active_rings,bool affine,
                 const std::vector<TemporalPrior>& priors={},double strength=1.) {
  Model model;constexpr int anchor=1;
  const bool fixed_anchor=priors.empty();
  std::array<bool,56> connected{};connected[anchor]=true;
  for(int j=0;j<56;++j)for(const auto& e:edges)if(e.ok && e.r<active_rings && (connected[e.r]||connected[e.s]))connected[e.r]=connected[e.s]=true;
  std::array<int,56> ai{},bi{};ai.fill(-1);bi.fill(-1);int nvar=0;
  for(int c=0;c<56;++c)if(connected[c]) {
    model.status[c]=fixed_anchor?(c==anchor?"fixed_energy_reference":"global_fit"):"global_fit_with_temporal_spectra";
    if(c!=anchor||!fixed_anchor){if(affine)ai[c]=nvar++;bi[c]=nvar++;}
  }
  if(nvar<2)return model;
  for(int iteration=0;iteration<30;++iteration) {
    TMatrixD normal(nvar,nvar);normal.Zero();TVectorD rhs(nvar);rhs.Zero();
    for(const auto& e:edges)if(e.ok && e.r<active_rings && connected[e.r] && connected[e.s]) {
      const double error=std::max(.05,e.scale/std::sqrt(std::max(1.,e.neff/2.)));
      for(double x:{e.xlo,e.xhi}) {
        const double y=e.a+e.b*x,residual=model.a[e.r]+model.b[e.r]*x-model.a[e.s]-model.b[e.s]*y;
        const double z=residual/std::max(.10,2.385*error),weight=1./(error*error*(iteration?1.+z*z:1.));
        std::vector<std::pair<int,double>> row;double constant=0;
        for(auto cq:{std::pair<int,double>{e.r,1.},std::pair<int,double>{e.s,-1.}}) {
          const int c=cq.first;const double sign=cq.second,energy=c==e.r?x:y;
          if(ai[c]>=0)row.emplace_back(ai[c],sign);
          if(bi[c]>=0)row.emplace_back(bi[c],sign*energy);else constant+=sign*energy;
        }
        for(const auto& iv:row) {
          rhs[iv.first]-=weight*iv.second*constant;
          for(const auto& jv:row)normal(iv.first,jv.first)+=weight*iv.second*jv.second;
        }
      }
    }
    // Independent channel spectra constrain the common scale and each response.
    // Their descriptive stability determines a weight, not a statistical prior.
    for(const auto& p:priors)if(connected[p.index]) {
      const int c=p.index;
      for(double x:{25.,75.}) {
        const double target=p.a+p.b*x,weight=strength/(p.error*p.error);
        std::vector<std::pair<int,double>> row;
        if(ai[c]>=0)row.emplace_back(ai[c],1.);
        if(bi[c]>=0)row.emplace_back(bi[c],x);
        for(const auto& iv:row) {
          rhs[iv.first]+=weight*iv.second*target;
          for(const auto& jv:row)normal(iv.first,jv.first)+=weight*iv.second*jv.second;
        }
      }
    }
    TDecompSVD solver(normal);Bool_t solved=false;TVectorD result=solver.Solve(rhs,solved);
    if(!solved)return model;
    const auto singular=solver.GetSig();double max=0,min=std::numeric_limits<double>::max();int rank=0;
    for(int i=0;i<singular.GetNrows();++i)max=std::max(max,singular[i]);
    for(int i=0;i<singular.GetNrows();++i)if(singular[i]>max*1e-12){min=std::min(min,singular[i]);++rank;}
    model.condition=max/min;model.rank=rank;model.iterations=iteration+1;
    if(rank!=nvar)return model;
    double change=0;
    for(int c=0;c<56;++c) {
      double a=ai[c]>=0?result[ai[c]]:0,b=bi[c]>=0?result[bi[c]]:1;
      change=std::max(change,std::abs(a-model.a[c])+50*std::abs(b-model.b[c]));model.a[c]=a;model.b[c]=b;
    }
    if(change<1e-6)break;
  }
  // Outer rings in high beam do not set the production sector calibration.
  // Fit each to all available calibrated sectors, preserving their diagnostics.
  for(int c=active_rings;c<24;++c) {
    double w=0,wx=0,wy=0,wxx=0,wxy=0;int count=0;
    for(const auto& e:edges)if(e.ok && e.r==c && connected[e.s]) {
      const double weight=1./std::pow(std::max(.05,e.scale/std::sqrt(std::max(1.,e.neff/2.))),2);
      for(double x:{e.xlo,e.xhi}) {
        double y=model.a[e.s]+model.b[e.s]*(e.a+e.b*x);
        w+=weight;wx+=weight*x;wy+=weight*y;wxx+=weight*x*x;wxy+=weight*x*y;++count;
      }
    }
    if(count>=4 && w*wxx-wx*wx>1e-9) {
      model.b[c]=affine?(w*wxy-wx*wy)/(w*wxx-wx*wx):wxy/wxx;
      model.a[c]=affine?(wy-model.b[c]*wx)/w:0.;model.status[c]="diagnostic_outer_ring_fit";
    }
  }
  model.ok=true;
  for(int c=0;c<56;++c)if(!(std::isfinite(model.a[c]) && std::isfinite(model.b[c]) && model.b[c]>.2 && model.b[c]<5. && std::abs(model.a[c])<50.)) {
    model.status[c]="unphysical_fit_reverted";model.a[c]=0;model.b[c]=1;model.ok=false;
  }
  return model;
}
json model_json(const Model& m,const json& channels) {
  json rows=json::array();
  for(const auto& c:channels) {
    const int i=c.at("index");const double olda=c.at("old_p0"),oldb=double(c.at("old_p1"))*double(c.at("old_att"));
    rows.push_back({{"index",i},{"channel",c.at("name")},{"status",m.status[i]},
      {"relative_a_keV",1000*m.a[i]},{"relative_b",m.b[i]},{"A_keV",1000*m.a[i]+m.b[i]*olda},{"B",m.b[i]*oldb}});
  }
  return {{"ok",m.ok},{"normal_matrix_condition",m.condition},{"rank",m.rank},{"iterations",m.iterations},{"coefficients",rows}};
}
// This constrained variant holds the total ADC-calibrated intercept fixed.
// Multiplying an already offset-corrected energy would also change its intercept,
// so the fitted slopes below act on E0 explicitly, in MeV.
Model fixed_offset_fit(const std::vector<Edge>& edges,int active,const json& channels,
                       const json& baseline,const json& temporal,double strength) {
  Model model;
  std::array<double,56> old_a{},old_b{},fixed_a{},slope{};
  std::array<bool,56> connected{};
  for(const auto& c:channels) {
    const int i=c.at("index");old_a[i]=double(c.at("old_p0"))/1000.;
    old_b[i]=double(c.at("old_p1"))*double(c.at("old_att"));
    require(old_b[i]>0.,"Nonpositive source calibration slope");
  }
  require(baseline.at("channels").size()==56,"Incomplete fixed-offset baseline");
  for(const auto& c:baseline.at("channels")) {
    const int i=c.at("index");require(i>=0&&i<56,"Invalid baseline index");
    fixed_a[i]=double(c.at("A_keV"))/1000.;slope[i]=c.at("B");
  }
  connected[1]=true;
  for(int j=0;j<56;++j)for(const auto& e:edges)
    if(e.ok&&e.r<active&&(connected[e.r]||connected[e.s]))connected[e.r]=connected[e.s]=true;
  std::array<int,56> index;index.fill(-1);int nvar=0;
  for(int c=0;c<56;++c)if(connected[c])index[c]=nvar++;
  if(nvar<2)return model;
  bool solved_all=false;
  for(int iteration=0;iteration<30;++iteration) {
    TMatrixD normal(nvar,nvar);normal.Zero();TVectorD rhs(nvar);rhs.Zero();
    auto add=[&](const std::vector<std::pair<int,double>>& row,double target,double weight) {
      for(const auto& i:row) {
        rhs[i.first]+=weight*i.second*target;
        for(const auto& j:row)normal(i.first,j.first)+=weight*i.second*j.second;
      }
    };
    for(const auto& e:edges)if(e.ok&&e.r<active&&connected[e.r]&&connected[e.s]) {
      const double error=std::max(.05,e.scale/std::sqrt(std::max(1.,e.neff/2.)));
      for(double x:{e.xlo,e.xhi}) {
        const double y=e.a+e.b*x,u=(x-old_a[e.r])/old_b[e.r],v=(y-old_a[e.s])/old_b[e.s];
        const double constant=fixed_a[e.r]-fixed_a[e.s];
        const double residual=constant+slope[e.r]*u-slope[e.s]*v;
        const double z=residual/std::max(.10,2.385*error);
        const double weight=1./(error*error*(iteration?1.+z*z:1.));
        add({{index[e.r],u},{index[e.s],-v}},-constant,weight);
      }
    }
    int priors=0;
    for(const auto& c:temporal.at("channels")) {
      const int i=c.at("index");
      if(!connected[i]||!c.at("training_supported").get<bool>())continue;
      ++priors;const double error=std::max(.5,50.*double(c.at("prior_log_scale")));
      for(double u:{25.,75.}) {
        const double target=double(c.at("direct_A_keV"))/1000.+double(c.at("direct_B"))*u-fixed_a[i];
        add({{index[i],u}},target,strength/(error*error));
      }
    }
    if(priors<3)return model;
    TDecompSVD solver(normal);Bool_t solved=false;TVectorD result=solver.Solve(rhs,solved);
    if(!solved)return model;
    const auto singular=solver.GetSig();double maximum=0,minimum=std::numeric_limits<double>::max();int rank=0;
    for(int i=0;i<singular.GetNrows();++i)maximum=std::max(maximum,singular[i]);
    for(int i=0;i<singular.GetNrows();++i)if(singular[i]>maximum*1e-12){minimum=std::min(minimum,singular[i]);++rank;}
    model.condition=maximum/minimum;model.rank=rank;model.iterations=iteration+1;
    if(rank!=nvar)return model;
    double change=0;
    for(int c=0;c<56;++c)if(connected[c]) {change=std::max(change,50*std::abs(slope[c]-result[index[c]]));slope[c]=result[index[c]];}
    solved_all=true;if(change<1e-6)break;
  }
  model.ok=solved_all;
  for(int c=0;c<56;++c) {
    model.b[c]=slope[c]/old_b[c];model.a[c]=fixed_a[c]-model.b[c]*old_a[c];
    model.status[c]=connected[c]?"fixed_total_offset_joint_slope":"baseline_disconnected";
    if(c>=active&&c<24)model.status[c]="diagnostic_outer_ring_baseline";
    if(!std::isfinite(model.a[c])||!std::isfinite(model.b[c])||model.b[c]<=.2||model.b[c]>=5.) {
      model.status[c]="unsupported_fixed_offset_fit";model.ok=false;
    }
  }
  return model;
}
json fit_temporal(const fs::path& edge_input,const fs::path& temporal_input,const fs::path& fixed_offset_input={}) {
  const json source=json::parse(read_text(edge_input)),temporal=json::parse(read_text(temporal_input));
  require(source.at("run")==temporal.at("run")&&source.at("period")==temporal.at("period"),"Temporal/edge input mismatch");
  std::vector<Edge> edges;
  for(const auto& row:source.at("edges").at("fit")) {
    Edge e;e.r=row.at("r");e.s=row.at("s");e.ok=row.at("ok");e.n=row.at("n");e.random=row.at("random");
    if(e.ok){e.a=row.at("a_MeV");e.b=row.at("b");e.scale=row.at("orthogonal_scale_MeV");e.neff=row.at("neff");e.xlo=row.at("xlo_MeV");e.xhi=row.at("xhi_MeV");}
    edges.push_back(e);
  }
  const int active=source.at("active_rings");std::vector<TemporalPrior> priors;
  for(const auto& c:temporal.at("channels"))if(c.at("training_supported").get<bool>()&&(int(c.at("index"))<active||int(c.at("index"))>=24)) {
    TemporalPrior p;p.index=c.at("index");
    const double olda=c.at("old_A_keV"),oldb=c.at("old_B");
    p.b=double(c.at("direct_B"))/oldb;p.a=(double(c.at("direct_A_keV"))-p.b*olda)/1000.;
    p.error=std::max(.5,50.*double(c.at("prior_log_scale")));priors.push_back(p);
  }
  json models=json::object(),prior_rows=json::array();
  for(const auto& p:priors)prior_rows.push_back({{"index",p.index},{"relative_a_MeV",p.a},{"relative_b",p.b},{"weight_scale_MeV",p.error}});
  if(!fixed_offset_input.empty()) {
    const auto document=json::parse(read_text(fixed_offset_input));
    const std::string run=source.at("run");const auto& baseline=document.at("runs").at(run);
    require(baseline.at("period")==source.at("period"),"Fixed-offset period mismatch");
    for(double strength:{.1,1.,10.,100.}) {
      const std::string name="fixed_offset_"+std::to_string(strength);
      models[name]=model_json(fixed_offset_fit(edges,active,source.at("sample_metadata").at("channels"),baseline,temporal,strength),source.at("sample_metadata").at("channels"));
      models[name]["temporal_weight_multiplier"]=strength;
      models[name]["fixed_total_offsets_keV"]=json::array();
      for(const auto& c:baseline.at("channels"))models[name]["fixed_total_offsets_keV"].push_back(c.at("A_keV"));
    }
  } else if(priors.size()>=3) {
    for(double strength:{.1,1.,10.,100.})for(bool affine:{false,true}) {
      std::string name=(affine?"affine_":"gain_")+std::to_string(strength);
      models[name]=model_json(global_fit(edges,active,affine,priors,strength),source.at("sample_metadata").at("channels"));
      models[name]["temporal_weight_multiplier"]=strength;
    }
  }
  json result={{"schema","sm150_s3_temporal_joint_fit_v1"},{"status","complete"},{"run",source.at("run")},{"period",source.at("period")},
    {"active_rings",active},{"edge_input",fs::absolute(edge_input).string()},{"edge_sha256",hash(read_text(edge_input))},
    {"temporal_input",fs::absolute(temporal_input).string()},{"temporal_sha256",hash(read_text(temporal_input))},
    {"source_pair_input",source.at("input")},{"priors",prior_rows},{"models",models},
    {"limitations","Only even-block data determine pixel lines and temporal estimates; finite stability weights are descriptive. Model selection and confirmation are separate."}};
  if(!fixed_offset_input.empty()) {
    result["fixed_offset_input"]=fs::absolute(fixed_offset_input).string();
    result["fixed_offset_sha256"]=hash(read_text(fixed_offset_input));
    result["fixed_offset_definition"]="Total A in E=A+B*E0 is held at the supplied baseline; slopes of connected production strips are solved jointly. Outer high-beam rings retain diagnostic baseline coefficients.";
  }
  return result;
}
}
int main(int argc,char** argv) {
  try {
    fs::path input,output,edge_input,temporal_input,fixed_offset_input;int minimum=80;
    for(int i=1;i<argc;++i) {
      const std::string arg=argv[i];auto value=[&](){require(++i<argc,"Missing value");return std::string(argv[i]);};
      if(arg=="--input")input=value();else if(arg=="--output")output=value();else if(arg=="--min-pairs")minimum=std::stoi(value());
      else if(arg=="--edge-input")edge_input=value();else if(arg=="--temporal-input")temporal_input=value();
      else if(arg=="--fixed-offset-baseline")fixed_offset_input=value();
      else if(arg=="--help"){std::cout<<"joint_calibrate --input pairs.root --output new.json [--min-pairs 80]\n"
        <<"joint_calibrate --edge-input fit.json --temporal-input conditional.json [--fixed-offset-baseline baseline.json] --output new.json\n";return 0;}
      else throw std::runtime_error("Unknown option: "+arg);
    }
    require(!output.empty() && !fs::exists(output) && !fs::exists(output.string()+".partial") && minimum>=20,"Invalid paths or existing output");
    if(!edge_input.empty()||!temporal_input.empty()) {
      require(input.empty()&&!edge_input.empty()&&!temporal_input.empty(),"Specify both temporal inputs, without --input");
      const auto begin=std::chrono::steady_clock::now();json result=fit_temporal(edge_input,temporal_input,fixed_offset_input);
      result["program_sha256"]=hash(read_text(fs::canonical("/proc/self/exe")));
      result["elapsed_seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
      fs::create_directories(fs::absolute(output).parent_path());fs::path partial=output.string()+".partial";
      std::ofstream stream(partial);stream<<result.dump(2)<<'\n';stream.close();require(bool(stream),"Cannot write temporal result");
      require(::link(partial.c_str(),output.c_str())==0&&::unlink(partial.c_str())==0,"Cannot publish temporal result");
      std::cout<<result.at("run")<<" temporal joint fit complete"<<std::endl;return 0;
    }
    require(!input.empty()&&fixed_offset_input.empty(),"Missing --input or fixed-offset mode without temporal inputs");
    const auto start=std::chrono::steady_clock::now();const Stat before(input);
    std::unique_ptr<TFile> file(TFile::Open(input.c_str(),"READ"));require(file && !file->IsZombie(),"Cannot read sample");
    auto* meta=dynamic_cast<TObjString*>(file->Get("metadata"));auto* tree=dynamic_cast<TTree*>(file->Get("pairs"));require(meta && tree,"Invalid sample");
    json metadata=json::parse(meta->GetString().Data());require(metadata.at("schema")=="sm150_s3_joint_sample_v1","Wrong schema");
    Int_t block=0,r=0,s=0,nr=0,ns=0,flags=0;Double_t x=0,y=0,dt=0;
    tree->SetBranchStatus("*",false);
    for(const auto& name:{"block","r","s","nr","ns","flags","r_old_keV","s_old_keV","dt_ns"})tree->SetBranchStatus(name,true);
    tree->SetBranchAddress("block",&block);tree->SetBranchAddress("r",&r);tree->SetBranchAddress("s",&s);tree->SetBranchAddress("nr",&nr);tree->SetBranchAddress("ns",&ns);
    tree->SetBranchAddress("flags",&flags);tree->SetBranchAddress("r_old_keV",&x);tree->SetBranchAddress("s_old_keV",&y);tree->SetBranchAddress("dt_ns",&dt);
    std::array<std::vector<Point>,768> points;std::array<int,768> sidebands{};uint64_t fit_pairs=0;
    for(Long64_t i=0;i<tree->GetEntries();++i) {
      require(tree->GetEntry(i)>0,"Failed pair read");require(r>=0 && r<24 && s>=24 && s<56,"Invalid strip index");
      if(block%2 || flags || nr!=1 || ns!=1 || x<15000 || x>=100000 || y<15000 || y>=100000)continue;
      const int pixel=r*32+s-24;
      if(dt>=-80 && dt<80){points[pixel].push_back({x/1000.,y/1000.,block});++fit_pairs;}
      else if((dt>=-880 && dt< -720)||(dt>=720 && dt<880))++sidebands[pixel];
    }
    json edge_sets=json::object(),models=json::object();const int active=metadata.at("period")=="low"?24:8;
    for(int subset:{-1,0,2}) {
      std::vector<Edge> edges;json erows=json::array();
      for(int p=0;p<768;++p) {
        auto e=fit_edge(points[p],p/32,24+p%32,sidebands[p],subset,minimum);
        edges.push_back(e);erows.push_back(edge_json(e));
      }
      const std::string name=subset<0?"fit":subset==0?"fit_early_alternating":"fit_late_alternating";
      edge_sets[name]=erows;
      for(bool affine:{false,true})models[name][affine?"affine":"gain"]=model_json(global_fit(edges,active,affine),metadata.at("channels"));
    }
    require(before.same(Stat(input)),"Sample changed while fitting");
    json result={{"schema","sm150_s3_joint_fit_v1"},{"status","complete"},{"run",metadata.at("run")},{"period",metadata.at("period")},{"input",fs::absolute(input).string()},
      {"sample_metadata",metadata},{"training_pairs",fit_pairs},{"active_rings",active},{"anchor_channel","ring02"},
      {"anchor_convention","Fix ring02 old corrected energy as gauge during relative fit; later common reference scaling must transform all strips together"},
      {"fit_selection","Even blocks, one hit per face above frozen 2.5MeV pre-FB threshold, clean flags, both old energies in [15,100)MeV, dt in [-80,80)ns"},
      {"random_selection","[-880,-720) and [720,880)ns; alpha=0.5; reject pixels where estimated true fraction <0.5"},
      {"minimum_pixel_pairs",minimum},{"minimum_energy_span_MeV",5.},{"edges",edge_sets},{"models",models},
      {"uncertainty","No full propagated coefficient uncertainty claimed; alternating training-block solutions record stability and gauge is exact by convention"},
      {"elapsed_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()},
      {"program_sha256",hash(read_text(fs::canonical("/proc/self/exe")))}};
    fs::create_directories(fs::absolute(output).parent_path());const fs::path partial=output.string()+".partial";
    std::ofstream out(partial);out<<result.dump(2)<<'\n';out.close();require(bool(out),"Cannot write result");
    require(::link(partial.c_str(),output.c_str())==0 && ::unlink(partial.c_str())==0,"Cannot publish result without overwriting");
    std::cout<<"complete run="<<metadata.at("run")<<" training_pairs="<<fit_pairs<<" output="<<output<<std::endl;return 0;
  }catch(const std::exception& e){std::cerr<<"joint S3 calibration: "<<e.what()<<'\n';return 1;}
}
