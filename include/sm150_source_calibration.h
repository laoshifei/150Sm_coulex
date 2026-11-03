#ifndef SM150_SOURCE_CALIBRATION_H
#define SM150_SOURCE_CALIBRATION_H
// Source-energy calibration utilities. No production parameters are installed.
#include "sm150_event_calibration.h"
#include <Fit/Fitter.h>
#include <Fit/FitResult.h>
#include <Math/Functor.h>
#include <TDecompSVD.h>
#include <TFile.h>
#include <TH1.h>
#include <TMatrixD.h>
#include <TObjString.h>
#include <TSpectrum.h>
#include <TVectorD.h>
#include <chrono>
#include <iostream>
#include <memory>
#include <numeric>
#include <set>

namespace sm150 { namespace sourcecal {
using diagnostic::json; using diagnostic::require; using diagnostic::read_text;
namespace fs=std::filesystem;
using Coeff=std::array<double,3>;
inline std::string sha_file(const fs::path& p) {
  std::ifstream in(p,std::ios::binary);require(bool(in),"Cannot read "+p.string());
  diagnostic::Digest digest;std::vector<char> b(1024*1024);
  while(in) {in.read(b.data(),b.size());if(in.gcount())digest.add(b.data(),size_t(in.gcount()));}
  require(in.eof(),"Hash read failed "+p.string());return digest.finish();
}
inline void write_json(const fs::path& p,const json& j) {
  std::ofstream out(p);out<<j.dump(2)<<'\n';out.close();require(bool(out),"Cannot write "+p.string());
}
inline double energy(double adc,const Coeff& p) {return p[0]+p[1]*adc+p[2]*adc*adc;}
inline double slope(double adc,const Coeff& p) {return p[1]+2*p[2]*adc;}
inline std::map<int,Coeff> read_seed(const fs::path& p,int cid) {
  std::istringstream input(read_text(p));std::string line;std::getline(input,line);
  std::map<int,Coeff> rows;std::set<int> seen;
  while(std::getline(input,line)) {
    if(line.empty())continue;
    std::istringstream row(line);int c,s,h;double a,b,d,chi;
    require(bool(row>>c>>s>>h>>a>>b>>d>>chi),"Invalid calibration row");
    require(c>=0&&c<=1&&s>=2&&s<=12&&h>=0&&h<16,"Invalid electronics address");
    const int key=c*10000+s*100+h;require(seen.insert(key).second,"Duplicate calibration address");
    require(std::isfinite(a)&&std::isfinite(b)&&std::isfinite(d),"Nonfinite seed");
    if(c!=cid || (a==0&&b==0&&d==0))continue;
    require(b>0 && slope(65535,{a,b,d})>0,"Nonmonotone seed");rows[key]={a,b,d};
  }
  require(!rows.empty(),"No physical seed channels");return rows;
}

struct Options {
  std::string command="source"; fs::path input,seed,output,config;
  int run=0;std::set<int> keys;bool legacy=false,help=false;
  Options(int argc,char** argv) {
    int i=1;if(i<argc && argv[i][0]!='-')command=argv[i++];
    for(;i<argc;++i) {
      std::string a=argv[i];
      if(a=="--help"||a=="-h"){help=true;continue;}
      if(a=="--legacy"){legacy=true;continue;}
      require(i+1<argc,"Missing value for "+a);std::string v=argv[++i];
      if(a=="--input")input=v;else if(a=="--seed")seed=v;else if(a=="--output")output=v;
      else if(a=="--config")config=v;else if(a=="--run")run=std::stoi(v);
      else if(a=="--keys") {
        std::istringstream parts(v);std::string item;while(std::getline(parts,item,','))keys.insert(std::stoi(item));
      } else throw std::runtime_error("Unknown option "+a);
    }
    if(help)return;
    require(!input.empty()&&!output.empty(),"Explicit --input and new --output are required");
    require(!fs::exists(output)&&!fs::is_symlink(output),"Refusing existing output directory "+output.string());
    input=fs::canonical(input);output=fs::absolute(output);
    if(command=="source"||command=="gain") {
      require(run>0&&!seed.empty(),"--run and --seed are required");seed=fs::canonical(seed);
    }
    if(config.empty())config=fs::canonical("/proc/self/exe").parent_path()/"source_config.json";
  }
};
inline void help(const char* kind) {
  std::cout<<kind<<"/calibrate source --run N --input HIST.root --seed TABLE.txt --output NEW_DIR [--keys K1,K2] [--config FILE] [--legacy]\n"
    <<"Modern input: chNNNNN/adc_clean; legacy input is explicit. Candidate files are never installed.\n";
}
struct Input {
  fs::path path; diagnostic::Stat before;std::unique_ptr<TFile> file;json metadata;
  explicit Input(const fs::path& p):path(p),before(p),file(TFile::Open(p.c_str(),"READ")) {
    require(file&&!file->IsZombie()&&!file->TestBit(TFile::kRecovered),"Invalid histogram ROOT "+p.string());
    if(auto* m=dynamic_cast<TObjString*>(file->Get("metadata")))metadata=json::parse(m->GetString().Data());
  }
  TH1* histogram(int key,const std::string& name,bool legacy) {
    std::ostringstream object;
    if(legacy)object<<"cid"<<key/10000<<"_sid"<<std::setw(2)<<std::setfill('0')<<key%10000/100
      <<"/h_cid"<<key/10000<<"_sid"<<std::setw(2)<<key%10000/100<<"_ch"<<std::setw(2)<<key%100;
    else object<<"ch"<<std::setw(5)<<std::setfill('0')<<key<<'/'<<name;
    auto* h=dynamic_cast<TH1*>(file->Get(object.str().c_str()));
    if(h)require(h->GetDimension()==1,"Expected one-dimensional source histogram");
    return h;
  }
  void verify() const {require(before.same(diagnostic::Stat(path)),"Input changed during calibration");}
};

// Deterministic affine matching in the energy coordinate defined by the seed.
// The seed is disclosed as prior matching information, not independent evidence.
inline json match_peaks(TH1& h,const Coeff& seed,const json& lines,const json& cfg) {
  const double emin=cfg.at("search_energy_keV")[0],emax=cfg.at("search_energy_keV")[1];
  auto copy=std::unique_ptr<TH1>(dynamic_cast<TH1*>(h.Clone("source_peak_search")));copy->SetDirectory(nullptr);
  int first=1,last=h.GetNbinsX();
  while(first<last&&energy(h.GetBinCenter(first),seed)<emin)++first;
  while(last>first&&energy(h.GetBinCenter(last),seed)>emax)--last;
  copy->GetXaxis()->SetRange(first,last);
  TSpectrum search(100);
  const int found=search.Search(copy.get(),cfg.value("search_sigma_bins",2.),"goff nodraw",cfg.at("search_threshold"));
  struct Peak {double adc,e,height;};std::vector<Peak> peaks;
  for(int i=0;i<found;++i) {
    const double adc=search.GetPositionX()[i]-.5,e=energy(adc,seed);
    if(e>=emin&&e<=emax)peaks.push_back({adc,e,search.GetPositionY()[i]});
  }
  const size_t limit=cfg.value("maximum_search_peaks",40);
  std::stable_sort(peaks.begin(),peaks.end(),[](const Peak&a,const Peak&b){return a.height>b.height;});
  if(peaks.size()>limit)peaks.resize(limit);
  std::sort(peaks.begin(),peaks.end(),[](const Peak&a,const Peak&b){return a.adc<b.adc;});
  std::vector<double> targets;for(const auto& l:lines)targets.push_back(l.at("energy_keV"));
  std::vector<int> best(targets.size(),-1);int best_count=-1;double best_score=1e100,besta=0,bestb=1;
  const double tol=cfg.at("matching_tolerance_keV"),offset_limit=cfg.at("matching_offset_limit_keV");
  const double bmin=cfg.at("matching_scale_range")[0],bmax=cfg.at("matching_scale_range")[1];
  auto evaluate=[&](double a,double b) {
    if(b<bmin||b>bmax||std::abs(a)>offset_limit)return;
    std::vector<int> selected(targets.size(),-1);int previous=-1,n=0;double error=0;
    for(size_t i=0;i<targets.size();++i) {
      double nearest=tol;int at=-1;
      for(size_t j=size_t(previous+1);j<peaks.size();++j) {
        const double d=std::abs(a+b*peaks[j].e-targets[i]);
        if(d<nearest){nearest=d;at=int(j);}
      }
      if(at>=0){selected[i]=at;previous=at;++n;error+=nearest*nearest/(tol*tol);}
    }
    const double score=error+.02*(a*a/(offset_limit*offset_limit)+(b-1)*(b-1)/(.2*.2));
    if(n>best_count||(n==best_count&&score<best_score)) {best_count=n;best_score=score;best=selected;besta=a;bestb=b;}
  };
  evaluate(0,1);
  for(size_t i=0;i<peaks.size();++i)for(size_t j=i+1;j<peaks.size();++j)
    for(size_t u=0;u<targets.size();++u)for(size_t v=u+1;v<targets.size();++v) {
      const double b=(targets[v]-targets[u])/(peaks[j].e-peaks[i].e);
      if(b<bmin||b>bmax)continue;
      evaluate(targets[u]-b*peaks[i].e,b);
    }
  json candidates=json::array(),matches=json::array();
  for(const auto& p:peaks)candidates.push_back({{"adc",p.adc},{"seed_energy_keV",p.e},{"height",p.height}});
  for(size_t i=0;i<targets.size();++i) {
    auto r=lines[i];r["matched"]=best[i]>=0;
    if(best[i]>=0){r["search_adc"]=peaks[best[i]].adc;r["search_seed_energy_keV"]=peaks[best[i]].e;}
    matches.push_back(r);
  }
  return {{"method","deterministic affine correction of seed coordinate; ordered one-to-one matching"},
          {"seed_used_for_matching",true},{"affine_offset_keV",besta},{"affine_scale",bestb},
          {"match_count",std::max(best_count,0)},{"search_candidates",candidates},{"lines",matches}};
}

inline double normal_cdf(double z) {return .5*std::erfc(-z/std::sqrt(2.));}
inline double log_normal_cdf(double z) {
  if(z>-10)return std::log(normal_cdf(z));
  const double t=1/(z*z);double series=1,term=1;
  for(int k=1;k<=10;++k){term*=-(2*k-1)*t;series+=term;}
  return -.5*z*z-std::log(-z)-.5*std::log(2*M_PI)+std::log(series);
}
inline double left_cdf(double x,double mu,double sigma,double tau) {
  const double z=(x-mu)/sigma;
  return std::clamp(normal_cdf(z)+std::exp((x-mu)/tau+.5*std::pow(sigma/tau,2)+log_normal_cdf(-z-sigma/tau)),0.,1.);
}
using Components=std::vector<std::pair<double,double>>;
inline std::array<double,2> bin_model(double l,double r,const double* p,bool tail,double lo,double hi,
                                    const Components& components={{0.,1.}},bool double_tail=false) {
  const double area=std::exp(p[0]),mu=p[1],sig=std::exp(p[2]);
  const double fraction=tail?p[5]:0.,tau=tail?std::exp(p[6]):1.;
  auto cdf=[&](double x){
    double sum=0;
    for(const auto& c:components) {
      if(double_tail)sum+=c.second*((1-fraction)*left_cdf(x,mu+c.first,sig,tau)+fraction*left_cdf(x,mu+c.first,sig,tau+std::exp(p[7])));
      else sum+=c.second*((1-fraction)*normal_cdf((x-mu-c.first)/sig)+fraction*left_cdf(x,mu+c.first,sig,tau));
    }
    return sum;
  };
  const double signal=area*std::max(0.,cdf(r)-cdf(l));
  const double bg=(r-l)*(std::exp(p[3])+(std::exp(p[4])-std::exp(p[3]))*((l+r)/2-lo)/(hi-lo));
  return {signal,bg};
}
inline json fit_peak(TH1& h,double guess,double halfwidth,bool tail,double physical_shift=-.5,
                     const Components& components={{0.,1.}},const json& policy=json::object(),
                     const json& start=json::object(),bool pure_tail=false,bool double_tail=false) {
  const double bin=h.GetBinWidth(1);
  int first=std::max(1,h.FindBin(guess-halfwidth-physical_shift));
  int last=std::min(h.GetNbinsX(),h.FindBin(guess+halfwidth-physical_shift));
  const int n=last-first+1,np=double_tail?8:tail?7:5,nfree=np-(pure_tail?1:0);
  require(n>nfree+3,"Insufficient bins for source peak fit");
  std::vector<double> edges,counts;
  double total=0,maximum=0;
  for(int i=first;i<=last;++i) {
    require(std::abs(h.GetBinWidth(i)-bin)<1e-8,"Nonuniform ADC histogram");
    const double y=h.GetBinContent(i);require(y>=0&&std::isfinite(y),"Invalid raw count bin");
    edges.push_back(h.GetBinLowEdge(i)+physical_shift);counts.push_back(y);total+=y;maximum=std::max(maximum,y);
  }
  edges.push_back(h.GetBinLowEdge(last+1)+physical_shift);require(total>=30,"Insufficient counts for source peak fit");
  std::vector<double> ends={counts[0],counts[1],counts[n-2],counts[n-1]};
  std::sort(ends.begin(),ends.end());const double background=std::max(.01,(ends[1]+ends[2])/(2*bin));
  const double area=std::max(20.,total-background*(edges.back()-edges.front()));
  std::vector<double> p={std::log(area),guess,std::log(std::max(bin*.7,halfwidth/8)),std::log(background),std::log(background)};
  std::vector<double> lower={std::log(.1),guess-halfwidth*.4,std::log(bin*.2),-18.,-18.};
  std::vector<double> upper={std::log(total*4+100),guess+halfwidth*.4,std::log(halfwidth*.45),std::log(maximum/bin*4+100),std::log(maximum/bin*4+100)};
  if(tail){p.push_back(start.value("fraction",.07));p.push_back(std::log(std::max(bin,halfwidth*start.value("tau_fraction",1./6))));lower.push_back(0);lower.push_back(std::log(bin*.15));upper.push_back(policy.value("tail_fraction_max",.6));upper.push_back(std::log(halfwidth));}
  if(double_tail){p.push_back(std::log(halfwidth*start.value("long_tau_increment_fraction",.3)));lower.push_back(std::log(bin*.15));upper.push_back(std::log(halfwidth));}
  p[1]+=start.value("center_shift_fraction",0.)*halfwidth;
  p[2]=std::log(std::max(bin*.7,halfwidth*start.value("sigma_fraction",.125)));
  for(int i=0;i<np;++i)p[i]=std::clamp(p[i],lower[i]+1e-6*(upper[i]-lower[i]),upper[i]-1e-6*(upper[i]-lower[i]));
  if(pure_tail)p[5]=1.;
  auto objective=[&](const double* q) {
    double sum=0;
    for(int i=0;i<n;++i) {
      const auto m=bin_model(edges[i],edges[i+1],q,tail,edges.front(),edges.back(),components,double_tail);const double v=std::max(m[0]+m[1],1e-200),y=counts[i];
      sum+=2*(v-y+(y>0?y*std::log(y/v):0.));
    }
    return std::max(0.,sum);
  };
  ROOT::Fit::Fitter fitter;
  fitter.Config().SetMinimizer("Minuit2","Migrad");
  fitter.Config().MinimizerOptions().SetPrintLevel(-1);
  fitter.Config().MinimizerOptions().SetMaxFunctionCalls(policy.value("maximum_function_calls",6000));
  fitter.Config().MinimizerOptions().SetTolerance(.001);
  fitter.Config().SetParamsSettings(np,p.data());
  for(int i=0;i<np;++i){fitter.Config().ParSettings(i).SetLimits(lower[i],upper[i]);fitter.Config().ParSettings(i).SetStepSize(i==1?bin*.1:.03);}
  if(pure_tail)fitter.Config().ParSettings(5).Fix();
  ROOT::Math::Functor functor(objective,np);
  const bool success=fitter.FitFCN(functor,p.data(),n,0);
  const auto& result=fitter.Result();
  for(int i=0;i<np;++i)p[i]=result.Parameter(i);
  json cov=json::array(),bounds=json::array();
  for(int i=0;i<np;++i) {
    json row=json::array();for(int j=0;j<np;++j)row.push_back(result.CovMatrix(i,j));cov.push_back(row);
    if(!(pure_tail&&i==5)&&std::min(p[i]-lower[i],upper[i]-p[i])<1e-4*(upper[i]-lower[i]))bounds.push_back(i);
  }
  json model=json::array(),signal=json::array(),bg=json::array(),residual=json::array();
  for(int i=0;i<n;++i) {
    const auto m=bin_model(edges[i],edges[i+1],p.data(),tail,edges.front(),edges.back(),components,double_tail);const double sum=std::max(m[0]+m[1],1e-200),y=counts[i];
    const double d=std::max(0.,2*(sum-y+(y>0?y*std::log(y/sum):0.)));
    model.push_back(sum);signal.push_back(m[0]);bg.push_back(m[1]);residual.push_back((y>=sum?1:-1)*std::sqrt(d));
  }
  const double dev=objective(p.data());bool shape_bound=false;
  for(const auto& i:bounds)if(i==1||i==2||i==0||
      ((i==5||i==6||i==7)&&p[i.get<int>()]>.5*(lower[i.get<int>()]+upper[i.get<int>()]))||(double_tail&&(i==5||i==6||i==7)))shape_bound=true;
  const double center_se=result.ParError(1),area_se=std::exp(p[0])*result.ParError(0);
  const bool usable=success&&result.IsValid()&&result.CovMatrixStatus()>=2&&!shape_bound&&
      std::isfinite(center_se)&&center_se>0&&center_se<std::exp(p[2])&&std::exp(p[0])>5*area_se;
  return {{"tail_model",double_tail?"double_left_tail":pure_tail?"pure_left_tail":tail?"left_tail":"gaussian"},{"success",success&&result.IsValid()},
    {"fit_status",result.Status()},{"covariance_status",result.CovMatrixStatus()},{"calls",result.NCalls()},
    {"center",p[1]},{"center_se",center_se},{"sigma",std::exp(p[2])},{"area",std::exp(p[0])},{"area_se",area_se},
    {"parameters",p},{"lower",lower},{"upper",upper},{"covariance",cov},{"bounds_indices",bounds},
    {"center_usable",usable},{"shape_bound",shape_bound},{"components_adc_offset_and_weight",components},
    {"deviance",dev},{"ndof",n-nfree},{"deviance_per_dof",dev/(n-nfree)},
    {"aicc",dev+2*nfree+2.*nfree*(nfree+1)/(n-nfree-1)},{"counts",counts},{"edges",edges},
    {"free_parameters",nfree},{"fixed_parameter_indices",pure_tail?json::array({5}):json::array()},
    {"initial_parameters",start},{"edm",result.Edm()},
    {"model",model},{"signal",signal},{"background",bg},{"residual",residual},
    {"error_scope","conditional Minuit covariance; no peak-shape or source-energy systematic uncertainty"}};
}
inline json choose_peak(TH1& h,double guess,double width,double shift=-.5,
                        const Components& components={{0.,1.}},const json& policy=json::object()) {
  json attempts=json::array();int selected=-1;double score=1e100;
  const int models=policy.value("include_double_left_tail",false)?4:policy.value("include_pure_left_tail",false)?3:2;
  for(int model=0;model<models;++model) {
    const bool tail=model>0,pure=model==2,dual=model==3;
    json starts=tail?policy.value("tail_starts",json::array({json::object()})):policy.value("gaussian_starts",json::array({json::object()}));
    json records=json::array(),best;double bestdev=1e100;bool bestusable=false;
    for(const auto& start:starts)try {
      auto fit=fit_peak(h,guess,width,tail,shift,components,policy,start,pure,dual);
      json record=fit;
      for(const auto* field:{"counts","edges","model","signal","background","residual"})record.erase(field);
      records.push_back(record);
      const bool usable=fit.at("center_usable");const double dev=fit.at("deviance");
      if(best.is_null()||(usable&&!bestusable)||(usable==bestusable&&dev<bestdev)) {
        best=std::move(fit);bestdev=dev;bestusable=usable;
      }
    } catch(const std::exception& e) {records.push_back({{"error",e.what()},{"initial_parameters",start},{"center_usable",false}});}
    if(best.is_null())best={{"tail_model",dual?"double_left_tail":pure?"pure_left_tail":tail?"left_tail":"gaussian"},{"center_usable",false}};
    best["start_records"]=records;
    if(best.value("center_usable",false)&&best.at("aicc").get<double>()<score){selected=int(attempts.size());score=best.at("aicc");}
    attempts.push_back(best);
  }
  if(selected<0)return {{"status","peak_fit_unreliable"},{"attempts",attempts}};
  json result=attempts[selected];result["status"]="conditional_peak_fit";result["attempts"]=attempts;result["selected_model"]=selected;
  return result;
}

inline json polynomial(const std::vector<json>& anchors,int degree) {
  const int n=int(anchors.size()),p=degree+1;require(n>=p,"Insufficient calibration anchors");
  double low=1e100,high=-1e100;for(const auto& a:anchors){const double x=a.at("center");low=std::min(low,x);high=std::max(high,x);}
  const double origin=(low+high)/2,scale=(high-low)/2;require(scale>0,"Degenerate calibration anchors");
  TMatrixD design(n,p);TVectorD target(n);
  for(int i=0;i<n;++i){const double z=(anchors[i].at("center").get<double>()-origin)/scale;for(int j=0;j<p;++j)design(i,j)=std::pow(z,j);target(i)=anchors[i].at("energy_keV");}
  Bool_t ok=false;TDecompSVD solver(design);const auto beta=solver.Solve(target,ok);require(ok,"Calibration SVD failed");
  Coeff coeff={beta[0]-beta[1]*origin/scale,beta[1]/scale,0};
  if(degree==2){coeff[2]=beta[2]/(scale*scale);coeff[1]-=2*origin*coeff[2];coeff[0]+=origin*origin*coeff[2];}
  TMatrixD transpose(TMatrixD::kTransposed,design),normal=transpose*design;
  TDecompSVD normal_solver(normal);auto inverse=normal_solver.Invert(ok);require(ok,"Calibration covariance solve failed");
  TMatrixD transform(3,p);transform.Zero();transform(0,0)=1;transform(0,1)=-origin/scale;transform(1,1)=1/scale;
  if(degree==2){transform(0,2)=origin*origin/(scale*scale);transform(1,2)=-2*origin/(scale*scale);transform(2,2)=1/(scale*scale);}
  TMatrixD response=transform*inverse*transpose,covariance(3,3);covariance.Zero();
  double ss=0,maximum=0,chi=0;json residual=json::array();
  for(int i=0;i<n;++i) {
    const double x=anchors[i].at("center"),y=anchors[i].at("energy_keV"),se=anchors[i].at("center_se");
    const double delta=energy(x,coeff)-y,var=std::pow(slope(x,coeff)*se,2);
    ss+=delta*delta;maximum=std::max(maximum,std::abs(delta));if(var>0)chi+=delta*delta/var;
    residual.push_back(delta);
    for(int a=0;a<3;++a)for(int b=0;b<3;++b)covariance(a,b)+=response(a,i)*var*response(b,i);
  }
  json cov=json::array();for(int a=0;a<3;++a){json row=json::array();for(int b=0;b<3;++b)row.push_back(covariance(a,b));cov.push_back(row);}
  return {{"coefficients",coeff},{"formula","E_keV=c0+c1*ADC+c2*ADC^2"},{"degree",degree},{"anchors",n},
    {"weights","equal anchors; covariance propagates conditional peak-center errors only"},{"covariance",cov},
    {"rms_keV",std::sqrt(ss/n)},{"max_residual_keV",maximum},{"conditional_center_chi2",chi},{"ndof",n-p},
    {"anchor_residuals_keV",residual},{"valid_adc",{low,high}},{"valid_energy_keV",{energy(low,coeff),energy(high,coeff)}},
    {"slope_at_zero",slope(0,coeff)},{"slope_at_65535",slope(65535,coeff)}};
}

inline json calibrate_channel(TH1& h,int key,const Coeff& seed,const json& lines,const json& cfg) {
  json result={{"key",key},{"seed_coefficients",seed},{"counts",h.Integral(0,h.GetNbinsX()+1)},
               {"independent_of_seed",false},{"formal_replacement",false}};
  auto matching=match_peaks(h,seed,lines,cfg);result["matching"]=matching;
  json all=json::array();std::vector<json> anchors;
  for(size_t i=0;i<matching["lines"].size();++i) {
    auto line=matching["lines"][i];line["used_for_calibration"]=false;
    if(!line.at("matched").get<bool>()){line["status"]="no_matching_peak";all.push_back(line);continue;}
    const double center=line.at("search_adc"),local_slope=slope(center,seed);
    double half=line.value("half_width_keV",cfg.value("peak_half_width_keV",8.))/local_slope;
    // Avoid fitting two source reference peaks as one local peak.
    for(size_t j=0;j<matching["lines"].size();++j)if(i!=j&&matching["lines"][j].at("matched").get<bool>())
      half=std::min(half,.44*std::abs(center-matching["lines"][j].at("search_adc").get<double>()));
    Components components={{0.,1.}};
    if(line.contains("components")) {
      components.clear();double total=0;
      const double scale=local_slope*matching.at("affine_scale").get<double>();
      for(const auto& c:line.at("components")) {
        const double w=c.at("probability_percent"),e=c.at("energy_keV");
        require(std::isfinite(w)&&w>0&&std::isfinite(e),"Invalid source multiplet");
        total+=w;components.emplace_back((e-line.at("energy_keV").get<double>())/scale,w);
      }
      require(total>0,"Empty source multiplet");for(auto& c:components)c.second/=total;
      line["fixed_component_scale_keV_per_adc"]=scale;
    }
    const auto fit=choose_peak(h,center,half,-.5,components,cfg.value("response_fit",json::object()));
    for(auto it=fit.begin();it!=fit.end();++it)line[it.key()]=it.value();
    if(line.contains("center")) {
      line["old_residual_keV"]=energy(line["center"],seed)-line.at("energy_keV").get<double>();
      const bool used=line.value("anchor",true)&&line.value("center_usable",false)&&
          line.at("deviance_per_dof").get<double>()<=cfg.at("maximum_peak_deviance_per_dof").get<double>();
      line["used_for_calibration"]=used;
      if(used)anchors.push_back(line);else line["anchor_rejection_reason"]=line.value("anchor",true)?"peak shape/fit quality":"predeclared diagnostic line";
    }
    all.push_back(line);
  }
  const int degree=cfg.at("degree"),minimum=cfg.at("minimum_anchors");
  result["candidate_eligible"]=false;
  if(int(anchors.size())>=degree+1) {
    auto curve=polynomial(anchors,degree);const auto p=curve.at("coefficients").get<Coeff>();
    bool eligible=int(anchors.size())>=minimum&&curve.at("max_residual_keV").get<double>()<=cfg.at("maximum_energy_residual_keV").get<double>()&&slope(0,p)>0&&slope(65535,p)>0;
    result["calibration"]=curve;result["candidate_eligible"]=eligible;
    result["status"]=eligible?(anchors.size()==size_t(degree+1)?"two_point_candidate_unvalidated":"candidate") : "calibration_review_required";
    for(auto& line:all)if(line.contains("center"))line["new_residual_keV"]=energy(line.at("center"),p)-line.at("energy_keV").get<double>();
  } else result["status"]="insufficient_reliable_anchors";
  result["peaks"]=all;return result;
}

inline int run_source(const Options& o,int cid,json (*channel_config)(int,const json&),
                     json (*channel_fit)(TH1&,int,const Coeff&,const json&,const json&)=calibrate_channel) {
  const auto start=std::chrono::steady_clock::now();const std::string config_text=read_text(o.config);
  const auto cfg=json::parse(config_text);const auto seeds=read_seed(o.seed,cid);Input input(o.input);
  if(!o.legacy){require(!input.metadata.is_null()&&input.metadata.value("schema","")=="sm150_raw2hist_v1","Expected raw2hist metadata");require(input.metadata.at("run").get<std::string>()==diagnostic::rn(o.run),"Source run does not match histogram metadata");}
  fs::create_directories(o.output);
  json result={{"schema","sm150_source_energy_candidate_v2"},{"run",diagnostic::rn(o.run)},
    {"input",o.input.string()},{"input_sha256",sha_file(o.input)},{"seed",o.seed.string()},{"seed_sha256",sha_file(o.seed)},
    {"config",cfg},{"config_sha256",diagnostic::hash(config_text)},{"input_metadata",input.metadata},
    {"binary_sha256",sha_file(fs::canonical("/proc/self/exe"))},{"formal_replacement",false},
    {"energy_convention","Source ADC calibration only; no in-beam hardware gain, run gain, attenuation or front/back terms"},
    {"channels",json::object()}};
  std::ofstream csv(o.output/"parameters.csv");csv<<"key,status,candidate_eligible,c0,c1,c2,anchors,rms_keV,max_residual_keV\n";
  std::ofstream text(o.output/"candidate_parameters.txt");text<<" cid sid ch p0 p1 p2 chi\n";text<<std::setprecision(17);csv<<std::setprecision(17);
  int selected=0,eligible=0;
  for(const auto& item:seeds) {
    const int key=item.first;if(!o.keys.empty()&&!o.keys.count(key))continue;++selected;
    json r;auto* h=input.histogram(key,"adc_clean",o.legacy);
    if(!h)r={{"key",key},{"status","missing_histogram"},{"candidate_eligible",false},{"seed_coefficients",item.second},{"peaks",json::array()}};
    else {
      try {auto settings=channel_config(key,cfg);r=channel_fit(*h,key,item.second,settings.at("lines"),settings);}
      catch(const std::exception& e){r={{"key",key},{"status","channel_exception"},{"error",e.what()},{"candidate_eligible",false},{"peaks",json::array()},{"seed_coefficients",item.second}};}
    }
    csv<<key<<','<<r.at("status").get<std::string>()<<','<<r.at("candidate_eligible").get<bool>();
    if(r.contains("calibration")) {
      const auto& c=r.at("calibration");const auto p=c.at("coefficients").get<Coeff>();
      csv<<','<<p[0]<<','<<p[1]<<','<<p[2]<<','<<c.at("anchors")<<','<<c.at("rms_keV")<<','<<c.at("max_residual_keV");
      if(r.at("candidate_eligible").get<bool>()) {++eligible;text<<key/10000<<' '<<key%10000/100<<' '<<key%100<<' '<<p[0]<<' '<<p[1]<<' '<<p[2]<<' '<<c.at("conditional_center_chi2")<<'\n';}
    } else csv<<",,,,,,";
    csv<<'\n';result["channels"][std::to_string(key)]=r;
    write_json(o.output/"calibration.json",result);
    std::cout<<diagnostic::rn(o.run)<<" channel "<<key<<": "<<r.at("status")<<std::endl;
  }
  require(selected>0 && (o.keys.empty()||selected==int(o.keys.size())),"Some requested keys are not physical seed channels");
  csv.close();text.close();require(bool(csv)&&bool(text),"Parameter output failed");input.verify();
  require(read_text(o.config)==config_text&&sha_file(o.seed)==result.at("seed_sha256"),"Calibration inputs changed");
  result["completion"]={{"status","computation_complete"},{"selected_channels",selected},{"candidate_channels",eligible},
      {"candidate_text_is_complete",selected==eligible},{"not_production_certified",true},
      {"elapsed_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()}};
  write_json(o.output/"calibration.json",result);write_json(o.output/"completion.json",result.at("completion"));
  return 0;
}
}}
#endif
