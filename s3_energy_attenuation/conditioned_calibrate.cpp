// Register signed, random-subtracted spectra at fixed ring-sector pixels.
// The cumulative-spectrum distance is descriptive, not a likelihood or chi-square.
#include "sm150_event_calibration.h"
#include <TFile.h>
#include <TH3I.h>
#include <TObjString.h>
#include <array>
#include <chrono>
#include <iostream>
#include <memory>
#include <numeric>
#include <unistd.h>
using namespace sm150::diagnostic;
namespace {
constexpr int bins=800,pixels=768,steps=121;
constexpr double width=250.,alpha=.5;
struct Cdf {
  std::array<double,bins+1> p{},r{};
  double at(const std::array<double,bins+1>& values,double energy)const {
    if(energy<=0)return 0;
    if(energy>=bins*width)return values.back();
    double v=energy/width;int i=int(v);return values[i]+(v-i)*(values[i+1]-values[i]);
  }
  double prompt(double lo,double hi)const{return at(p,hi)-at(p,lo);}
  double random(double lo,double hi)const{return at(r,hi)-at(r,lo);}
  double net(double lo,double hi)const{return prompt(lo,hi)-alpha*random(lo,hi);}
};
struct Data {
  json meta;std::vector<Cdf> counts;
  Data():counts(2*pixels*4){}
  const Cdf& get(int face,int pixel,int split)const{return counts[(face*pixels+pixel)*4+split];}
  Cdf& get(int face,int pixel,int split){return counts[(face*pixels+pixel)*4+split];}
};
Data load(const fs::path& path) {
  Data d;std::unique_ptr<TFile> f(TFile::Open(path.c_str(),"READ"));require(f&&!f->IsZombie(),"Cannot read spectra");
  auto* m=dynamic_cast<TObjString*>(f->Get("metadata"));require(m,"Missing metadata");d.meta=json::parse(m->GetString().Data());
  require(d.meta.at("schema")=="sm150_s3_conditioned_spectra_v1","Wrong schema");
  for(int face=0;face<2;++face) {
    std::string name=face?"sector":"ring";
    auto* p=dynamic_cast<TH3I*>(f->Get((name+"_prompt").c_str()));auto* r=dynamic_cast<TH3I*>(f->Get((name+"_random").c_str()));
    require(p&&r&&p->GetNbinsX()==bins&&p->GetNbinsY()==pixels&&p->GetNbinsZ()==4,"Wrong histogram");
    for(int pixel=0;pixel<pixels;++pixel)for(int split=0;split<4;++split) {
      auto& c=d.get(face,pixel,split);
      for(int i=0;i<bins;++i){c.p[i+1]=c.p[i]+p->GetBinContent(i+1,pixel+1,split+1);c.r[i+1]=c.r[i]+r->GetBinContent(i+1,pixel+1,split+1);}
    }
  }
  return d;
}
Cdf combine(const Data& d,int face,int pixel,int split) {
  if(split>=0)return d.get(face,pixel,split);
  Cdf c;const auto& a=d.get(face,pixel,0);const auto& b=d.get(face,pixel,2);
  for(int i=0;i<=bins;++i){c.p[i]=a.p[i]+b.p[i];c.r[i]=a.r[i]+b.r[i];}return c;
}
struct Comparison {
  Cdf reference,current;double ar=0,br=1,ac=0,bc=1,lo=20000,hi=90000;
  std::array<double,3> counts(bool ref,double q=1.)const {
    const auto& c=ref?reference:current;double a=ref?ar:q*ac,b=ref?br:q*bc;
    double x=(lo-a)/b,y=(hi-a)/b;
    if(!(b>0&&x>=0&&y<=200000))return {0,0,0};
    double p=c.prompt(x,y),r=c.random(x,y);return {p-alpha*r,p,r};
  }
  double distance_map(double a,double b)const {
    double rx=(lo-ar)/br,ry=(hi-ar)/br,cx=(lo-a)/b,cy=(hi-a)/b;
    if(!(br>0&&b>0&&rx>=0&&ry<=200000&&cx>=0&&cy<=200000))return 1.;
    double nr=reference.net(rx,ry),nc=current.net(cx,cy);
    if(nr<15||nc<15)return 1.;
    double total=0;int n=0;
    for(double e=lo+1000;e<hi;e+=1000) {
      double r=reference.net(rx,(e-ar)/br)/nr,c=current.net(cx,(e-a)/b)/nc;
      total+=(c-r)*(c-r);++n;
    }
    return std::sqrt(total/std::max(1,n));
  }
  double operator()(double q)const{return distance_map(q*ac,q*bc);}
  bool usable()const {
    auto r=counts(true),c=counts(false);
    return r[0]>=40&&c[0]>=40&&r[0]>=.5*r[1]&&c[0]>=.5*c[1];
  }
};
struct Estimate{double q=1,value=1,low=.5,high=2;bool boundary=false;std::vector<double> grid,profile;};
template<class F> Estimate estimate(F fn) {
  Estimate v;int best=0;
  for(int i=0;i<steps;++i) {
    double q=std::exp(std::log(.5)+i*std::log(4.)/(steps-1));v.grid.push_back(q);v.profile.push_back(fn(q));
    if(v.profile.back()<v.profile[best])best=i;
  }
  v.boundary=best==0||best==steps-1;
  double lo=v.grid[std::max(0,best-1)],hi=v.grid[std::min(steps-1,best+1)],phi=(std::sqrt(5.)-1)/2;
  double c=hi-phi*(hi-lo),d=lo+phi*(hi-lo),fc=fn(c),fd=fn(d);
  for(int i=0;i<28;++i) {
    if(fc<fd){hi=d;d=c;fd=fc;c=hi-phi*(hi-lo);fc=fn(c);}
    else{lo=c;c=d;fc=fd;d=lo+phi*(hi-lo);fd=fn(d);}
  }
  v.q=.5*(lo+hi);v.value=fn(v.q);std::vector<double> support{v.q};
  for(int i=0;i<steps;++i)if(v.profile[i]<=v.value+.005)support.push_back(v.grid[i]);
  v.low=*std::min_element(support.begin(),support.end());v.high=*std::max_element(support.begin(),support.end());
  double outer_low=.5,outer_high=2.;
  for(double q:v.grid){if(q<v.low)outer_low=q;if(q>v.high){outer_high=q;break;}}
  double inside=v.low,outside=outer_low;
  for(int i=0;i<24;++i){double mid=.5*(inside+outside);if(fn(mid)<=v.value+.005)inside=mid;else outside=mid;}v.low=inside;
  inside=v.high;outside=outer_high;
  for(int i=0;i<24;++i){double mid=.5*(inside+outside);if(fn(mid)<=v.value+.005)inside=mid;else outside=mid;}v.high=inside;
  return v;
}
double group_distance(const std::vector<Comparison>& groups,double q) {
  if(groups.empty())return 1.;
  double total=0;for(const auto& g:groups){double d=g(q);total+=d*d;}return std::sqrt(total/groups.size());
}
double median(std::vector<double> v) {
  if(v.empty())return 0.;
  std::sort(v.begin(),v.end());size_t n=v.size();return n%2?v[n/2]:.5*(v[n/2-1]+v[n/2]);
}
}
int main(int argc,char** argv) {
  try {
    fs::path input,reference,baseline,output;
    for(int i=1;i<argc;++i) {
      std::string arg=argv[i];auto value=[&](){require(++i<argc,"Missing argument");return std::string(argv[i]);};
      if(arg=="--input")input=value();else if(arg=="--reference")reference=value();else if(arg=="--baseline")baseline=value();else if(arg=="--output")output=value();
      else throw std::runtime_error("Unknown option "+arg);
    }
    require(!input.empty()&&!reference.empty()&&!baseline.empty()&&!output.empty()&&!fs::exists(output)&&!fs::exists(output.string()+".partial"),"Missing/existing paths");
    auto start=std::chrono::steady_clock::now();Data current=load(input),ref=load(reference);
    std::string run=current.meta.at("run"),ref_run=ref.meta.at("run"),period=current.meta.at("period");require(period==ref.meta.at("period"),"Period mismatch");
    int active=period=="low"?24:8;
    json package=json::parse(read_text(baseline));auto base=package.at("runs").at(run).at("channels"),rb=package.at("runs").at(ref_run).at("channels");
    std::array<double,56> a{},b{},ra{},rbb{},oa{},ob{};
    for(const auto& c:base){int i=c.at("index");a[i]=c.at("A_keV");b[i]=c.at("B");oa[i]=c.at("old_A_keV");ob[i]=c.at("old_B");}
    for(const auto& c:rb){int i=c.at("index");ra[i]=c.at("A_keV");rbb[i]=c.at("B");}
    json rows=json::array(),pixel_rows=json::array();
    for(int channel=0;channel<56;++channel) {
      int face=channel<24?0:1;std::vector<int> eligible;std::array<std::vector<Comparison>,5> sets;
      // The last set holds the combined even blocks; indices0..3 are actual folds.
      auto make=[&](int pixel,int split){Comparison c;c.reference=combine(ref,face,pixel,split);c.current=combine(current,face,pixel,split);c.ar=ra[channel];c.br=rbb[channel];c.ac=a[channel];c.bc=b[channel];return c;};
      std::vector<int> all_pixels;
      if(face==0)for(int s=0;s<32;++s)all_pixels.push_back(channel*32+s);
      else for(int r=0;r<active;++r)all_pixels.push_back(r*32+channel-24);
      double total_current=0,total_ref=0;
      for(int pixel:all_pixels) {
        auto train=make(pixel,-1);auto cn=train.counts(false),rn=train.counts(true);
        if(train.usable()) {
          eligible.push_back(pixel);total_current+=cn[0];total_ref+=rn[0];sets[4].push_back(train);
          for(int split=0;split<4;++split)sets[split].push_back(make(pixel,split));
        }
      }
      auto fit=estimate([&](double q){return group_distance(sets[4],q);});
      auto early=estimate([&](double q){return group_distance(sets[0],q);});
      auto late=estimate([&](double q){return group_distance(sets[2],q);});
      auto alternate=sets[4];for(auto& g:alternate){g.lo=30000;g.hi=80000;}
      auto alt=estimate([&](double q){return group_distance(alternate,q);});
      std::vector<std::string> training_flags;
      if(eligible.size()<3)training_flags.push_back("fewer_than_three_usable_pixels");
      if(total_current<300||total_ref<300)training_flags.push_back("net_training_count_below_300");
      if(fit.boundary)training_flags.push_back("search_boundary");
      double split=std::abs(std::log(early.q/late.q)),window=std::abs(std::log(alt.q/fit.q)),support=std::log(fit.high/fit.low);
      if(split>.05)training_flags.push_back("alternating_fit_difference_above5pct");
      if(window>.05)training_flags.push_back("energy_window_difference_above5pct");
      if(support>.20)training_flags.push_back("broad_descriptive_profile");
      std::vector<double> penalties,q_pixels;int worse=0;
      for(int pixel:all_pixels) {
        auto c=make(pixel,-1);bool valid=c.usable();json pr={{"channel_index",channel},{"ring",pixel/32+1},{"sector",pixel%32+1},{"usable",valid},
          {"training_counts",c.counts(false)},{"reference_counts",c.counts(true)}};
        if(valid) {
          auto single=estimate([&](double q){return c(q);});double shared=c(fit.q);penalties.push_back(shared-single.value);q_pixels.push_back(single.q);
          if(shared>c(1.)+.015)++worse;
          auto sel=make(pixel,1),check=make(pixel,3);
          pr.update({{"q",single.q},{"boundary",single.boundary},{"fit_D_base",c(1.)},{"fit_D_individual",single.value},{"fit_D_shared",shared},
            {"selection_D_base",sel(1.)},{"selection_D_shared",sel(fit.q)},{"selection_D_individual",sel(single.q)},
            {"confirmation_D_base",check(1.)},{"confirmation_D_shared",check(fit.q)},{"confirmation_D_individual",check(single.q)},
            {"selection_counts",sel.counts(false,fit.q)},{"confirmation_counts",check.counts(false,fit.q)}});
        }
        pixel_rows.push_back(pr);
      }
      double incompat=median(penalties);
      if(incompat>.025)training_flags.push_back("shared_gain_loses_pixel_shape");
      if(worse>int(.25*eligible.size()))training_flags.push_back("shared_gain_worsens_many_pixels");
      std::vector<double> variations;for(double q:q_pixels)variations.push_back(std::abs(std::log(q/fit.q)));
      json metrics=json::object();
      for(int fold:{1,3,4}) {
        double old=0;for(const auto& g:sets[fold])old+=std::pow(g.distance_map(oa[channel],ob[channel]),2);
        metrics[fold==1?"selection":fold==3?"confirmation":"fit"]={{"D_formal",sets[fold].empty()?1.:std::sqrt(old/sets[fold].size())},
          {"D_0228",group_distance(sets[fold],1.)},{"D_direct",group_distance(sets[fold],fit.q)}};
      }
      std::vector<std::string> selection_flags,confirmation_flags;
      double ds=metrics["selection"]["D_0228"],dn=metrics["selection"]["D_direct"];
      if(dn>.9*ds||ds-dn<.002)selection_flags.push_back("no_clear_selection_gain");
      if(dn>.12)selection_flags.push_back("large_selection_shape_residual");
      if(double(metrics["confirmation"]["D_direct"])>double(metrics["confirmation"]["D_0228"])+.002)confirmation_flags.push_back("confirmation_shape_not_improved");
      bool fixed=run==ref_run;if(fixed){fit.q=1.;training_flags.clear();selection_flags.clear();confirmation_flags.clear();}
      double uncertainty=std::max({.01,.5*split,.5*window,.5*support});
      rows.push_back({{"index",channel},{"channel",base[channel].at("channel")},{"production_used",channel>=24||channel<active},
        {"q",fit.q},{"q_fit_fold0",early.q},{"q_fit_fold2",late.q},{"q_alternate_window",alt.q},{"q_profile_support",{fit.low,fit.high}},
        {"q_profile",{{"grid",fit.grid},{"distance",fit.profile}}},{"eligible_pixels",eligible},{"total_training_net",total_current},{"reference_training_net",total_ref},
        {"pixel_log_spread_median",median(variations)},{"shared_pixel_penalty_median",incompat},{"split_log_difference",split},{"window_log_difference",window},
        {"training_flags",training_flags},{"selection_flags",selection_flags},{"confirmation_flags",confirmation_flags},
        {"training_supported",training_flags.empty()&&!fixed},{"fixed_reference",fixed},{"metrics",metrics},
        {"base_A_keV",a[channel]},{"base_B",b[channel]},{"old_A_keV",oa[channel]},{"old_B",ob[channel]},
        {"reference_A_keV",ra[channel]},{"reference_B",rbb[channel]},
        {"direct_A_keV",fit.q*a[channel]},{"direct_B",fit.q*b[channel]},
        {"prior_log_scale",uncertainty},{"uncertainty_status","Descriptive stability scale, not a confidence interval"}});
    }
    json result={{"schema","sm150_s3_conditioned_fit_v1"},{"status","complete"},{"run",run},{"period",period},{"reference_run",ref_run},
      {"input",fs::absolute(input).string()},{"reference",fs::absolute(reference).string()},{"input_sha256",hash(read_text(input))},
      {"reference_sha256",hash(read_text(reference))},{"baseline_sha256",hash(read_text(baseline))},{"channels",rows},{"pixels",pixel_rows},
      {"distance","RMS of difference of separately normalized signed net cumulative spectra; P-0.5R, no clipping; equal eligible-pixel weights,1MeV grid20-90MeV"},
      {"program_sha256",hash(read_text(fs::canonical("/proc/self/exe")))},{"elapsed_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()}};
    fs::create_directories(fs::absolute(output).parent_path());fs::path partial=output.string()+".partial";
    std::ofstream stream(partial);stream<<result.dump(2)<<'\n';stream.close();require(bool(stream),"Write failed");
    require(::link(partial.c_str(),output.c_str())==0&&::unlink(partial.c_str())==0,"Publication failed");
    std::cout<<run<<" conditioned fit complete"<<std::endl;return 0;
  }catch(const std::exception& e){std::cerr<<"conditioned fit: "<<e.what()<<'\n';return 1;}
}
