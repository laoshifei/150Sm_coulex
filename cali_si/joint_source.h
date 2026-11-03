#ifndef SM150_JOINT_ALPHA_SOURCE_H
#define SM150_JOINT_ALPHA_SOURCE_H
#include "../include/sm150_source_calibration.h"

namespace sm150 { namespace sourcecal {

// Fit all source features together: the tail of an upper feature contributes
// to every lower-energy bin. Local straight backgrounds cannot replace it.
inline json joint_alpha_channel(TH1& h,int key,const Coeff& seed,const json& lines,const json& cfg) {
  if(!cfg.value("joint_alpha_fit",false))return calibrate_channel(h,key,seed,lines,cfg);
  auto matching=match_peaks(h,seed,lines,cfg);
  std::vector<json> refs;std::vector<Components> components;
  std::vector<double> guesses,halfwidths,starts,lower,upper;
  const double bin=h.GetBinWidth(1);
  for(const auto& line:matching.at("lines"))if(line.at("matched").get<bool>()) {
    refs.push_back(line);guesses.push_back(line.at("search_adc"));
    halfwidths.push_back(cfg.at("peak_half_width_keV").get<double>()/slope(guesses.back(),seed));
    Components comp;double total=0.;
    const double scale=slope(guesses.back(),seed)*matching.at("affine_scale").get<double>();
    for(const auto& c:line.at("components")) {
      const double weight=c.at("probability_percent");total+=weight;
      comp.emplace_back((c.at("energy_keV").get<double>()-line.at("energy_keV").get<double>())/scale,weight);
    }
    for(auto& c:comp)c.second/=total;
    components.push_back(comp);
  }
  json result={{"key",key},{"seed_coefficients",seed},{"counts",h.Integral(0,h.GetNbinsX()+1)},
    {"independent_of_seed",false},{"formal_replacement",false},{"matching",matching},
    {"candidate_eligible",false},{"status","insufficient_reliable_anchors"},{"peaks",json::array()}};
  const bool dual=cfg.value("joint_response",std::string("gaussian_left_mixture"))=="double_left_tail";
  const int stride=dual?6:5;
  const int count=int(refs.size());
  if(count<2){result["status"]="insufficient_matched_source_features";return result;}
  const int first=std::max(1,h.FindBin(guesses.front()-halfwidths.front()+.5));
  const int last=std::min(h.GetNbinsX(),h.FindBin(guesses.back()+halfwidths.back()+.5));
  const int n=last-first+1,np=stride*count+2;
  require(n>np+3,"Insufficient bins for joint source fit");
  std::vector<double> edges,counts;
  double maximum=0.,total=0.;
  for(int b=first;b<=last;++b) {
    const double y=h.GetBinContent(b);require(y>=0&&std::isfinite(y),"Invalid source count");
    edges.push_back(h.GetBinLowEdge(b)-.5);counts.push_back(y);maximum=std::max(maximum,y);total+=y;
  }
  edges.push_back(h.GetBinLowEdge(last+1)-.5);
  require(total>=100,"Insufficient source counts");
  for(int k=0;k<count;++k) {
    double area=0.;for(int i=0;i<n;++i)if(std::abs((edges[i]+edges[i+1])/2-guesses[k])<halfwidths[k])area+=counts[i];
    starts.insert(starts.end(),{std::log(std::max(30.,area)),guesses[k]+.05*halfwidths[k],std::log(std::max(.8*bin,halfwidths[k]*.08)),.35,std::log(halfwidths[k]*.2)});
    lower.insert(lower.end(),{std::log(.1),guesses[k]-.45*halfwidths[k],std::log(.2*bin),0.,std::log(.15*bin)});
    upper.insert(upper.end(),{std::log(total*4+100),guesses[k]+.45*halfwidths[k],std::log(.45*halfwidths[k]),1.,std::log(halfwidths[k]*1.5)});
    if(dual){starts[stride*k+4]=std::log(halfwidths[k]*.04);starts.push_back(std::log(halfwidths[k]*.35));lower.push_back(std::log(.15*bin));upper.push_back(std::log(halfwidths[k]*1.5));}
  }
  const double bg=std::max(.01,(counts.front()+counts.back())/(2*bin));
  starts.insert(starts.end(),{std::log(bg),std::log(bg)});
  lower.insert(lower.end(),{-18.,-18.});upper.insert(upper.end(),{std::log(maximum/bin*4+100),std::log(maximum/bin*4+100)});
  auto feature_cdf=[&](double x,const double* q,int k) {
    const double mu=q[stride*k+1],sigma=std::exp(q[stride*k+2]),fraction=q[stride*k+3],tau=std::exp(q[stride*k+4]);
    double value=0.;for(const auto& c:components[k]) {
      if(dual)value+=c.second*((1-fraction)*left_cdf(x,mu+c.first,sigma,tau)+fraction*left_cdf(x,mu+c.first,sigma,std::exp(q[stride*k+5])));
      else value+=c.second*((1-fraction)*normal_cdf((x-mu-c.first)/sigma)+fraction*left_cdf(x,mu+c.first,sigma,tau));
    }
    return value;
  };
  auto models=[&](const double* q,std::vector<std::vector<double>>* separate=nullptr) {
    std::vector<double> model(n);
    for(int i=0;i<n;++i)model[i]=(edges[i+1]-edges[i])*(std::exp(q[stride*count])+
      (std::exp(q[stride*count+1])-std::exp(q[stride*count]))*((edges[i]+edges[i+1])/2-edges.front())/(edges.back()-edges.front()));
    if(separate){separate->clear();separate->push_back(model);}
    for(int k=0;k<count;++k) {
      std::vector<double> part(n);double previous=feature_cdf(edges.front(),q,k);const double area=std::exp(q[stride*k]);
      for(int i=0;i<n;++i){const double current=feature_cdf(edges[i+1],q,k);part[i]=area*std::max(0.,current-previous);model[i]+=part[i];previous=current;}
      if(separate)separate->push_back(std::move(part));
    }
    return model;
  };
  auto objective=[&](const double* q) {
    const auto model=models(q);double dev=0.;
    for(int i=0;i<n;++i){const double m=std::max(model[i],1e-200),y=counts[i];dev+=2*(m-y+(y>0?y*std::log(y/m):0.));}
    return std::max(0.,dev);
  };
  auto optimize=[&](std::vector<double> initial,const std::set<int>& fixed) {
    ROOT::Fit::Fitter fitter;fitter.Config().SetMinimizer("Minuit2","Migrad");
    fitter.Config().MinimizerOptions().SetPrintLevel(-1);
    fitter.Config().MinimizerOptions().SetMaxFunctionCalls(cfg.value("joint_maximum_calls",18000));
    fitter.Config().MinimizerOptions().SetTolerance(.001);fitter.Config().SetParamsSettings(np,initial.data());
    for(int i=0;i<np;++i){fitter.Config().ParSettings(i).SetLimits(lower[i],upper[i]);fitter.Config().ParSettings(i).SetStepSize(i<stride*count&&i%stride==1?.1*bin:.03);if(fixed.count(i))fitter.Config().ParSettings(i).Fix();}
    ROOT::Math::Functor functor(objective,np);const bool ok=fitter.FitFCN(functor,initial.data(),n,0);const auto& fit=fitter.Result();
    json covariance=json::array(),errors=json::array(),bounds=json::array();
    for(int i=0;i<np;++i){initial[i]=fit.Parameter(i);errors.push_back(fit.ParError(i));json row=json::array();for(int j=0;j<np;++j)row.push_back(fit.CovMatrix(i,j));covariance.push_back(row);
      if(!fixed.count(i)&&std::min(initial[i]-lower[i],upper[i]-initial[i])<1e-4*(upper[i]-lower[i]))bounds.push_back(i);}
    return json{{"success",ok&&fit.IsValid()},{"fit_status",fit.Status()},{"covariance_status",fit.CovMatrixStatus()},
      {"parameters",initial},{"errors",errors},{"covariance",covariance},{"deviance",objective(initial.data())},
      {"edm",fit.Edm()},{"calls",fit.NCalls()},{"fixed_indices",fixed},{"bounds_indices",bounds},
      {"ndof",n-np+int(fixed.size())}};
  };
  json attempts=json::array(),best;double bestdev=1e100;bool bestok=false;
  const std::vector<double> fractions={.2,.65,.95};
  for(double fraction:fractions) {
    auto initial=starts;for(int k=0;k<count;++k){initial[stride*k+3]=fraction;initial[stride*k+1]=guesses[k]+fraction*.1*halfwidths[k];}
    auto fit=optimize(initial,{});fit["start_fraction"]=fraction;attempts.push_back(fit);
    const bool ok=fit.at("success")&&fit.at("covariance_status").get<int>()>=2;const double dev=fit.at("deviance");
    if(best.is_null()||(ok&&!bestok)||(ok==bestok&&dev<bestdev)){best=fit;bestdev=dev;bestok=ok;}
  }
  // Nested endpoint models remove unidentifiable parameters, rather than
  // treating a numerically saturated fraction as a measured interior value.
  auto p=best.at("parameters").get<std::vector<double>>();std::set<int> fixed;
  for(int k=0;k<count;++k) {
    if(dual) {
      if(p[stride*k+4]>p[stride*k+5]){std::swap(p[stride*k+4],p[stride*k+5]);p[stride*k+3]=1.-p[stride*k+3];}
      if(p[stride*k+3]>.998){p[stride*k+4]=p[stride*k+5];p[stride*k+3]=0.;fixed.insert(stride*k+3);fixed.insert(stride*k+5);}
      else if(p[stride*k+3]<.002){p[stride*k+3]=0.;fixed.insert(stride*k+3);fixed.insert(stride*k+5);}
    } else {
      if(p[stride*k+3]>.998){p[stride*k+3]=1.;fixed.insert(stride*k+3);}
      else if(p[stride*k+3]<.002||p[stride*k+4]<lower[stride*k+4]+.01){p[stride*k+3]=0.;fixed.insert(stride*k+3);fixed.insert(stride*k+4);}
    }
  }
  if(dual||!fixed.empty()) {
    auto reduced=optimize(p,fixed);reduced["endpoint_reduction"]=true;attempts.push_back(reduced);
    if(reduced.at("success").get<bool>()&&reduced.at("covariance_status").get<int>()>=2&&reduced.at("deviance").get<double>()<=bestdev+2*fixed.size())best=reduced;
  }
  p=best.at("parameters").get<std::vector<double>>();const auto errors=best.at("errors").get<std::vector<double>>();
  fixed=best.at("fixed_indices").get<std::set<int>>();std::set<int> bound=best.at("bounds_indices").get<std::set<int>>();
  std::vector<std::vector<double>> parts;const auto model=models(p.data(),&parts);
  result["joint_fit"]=best;result["joint_fit"]["attempts"]=attempts;result["joint_fit"]["edges"]=edges;
  result["joint_fit"]["counts"]=counts;result["joint_fit"]["model"]=model;result["joint_fit"]["components"]=parts;
  result["joint_fit"]["lower"]=lower;result["joint_fit"]["upper"]=upper;
  result["joint_fit"]["response"]=dual?"Sum of all alpha multiplets with per-feature two-exponential left response convolved with a Gaussian; shared linear background":"Sum of all alpha multiplets with per-feature Gaussian/left-exGaussian response and shared linear background";
  result["joint_fit"]["parameter_stride"]=stride;
  std::vector<json> anchors;
  for(int k=0;k<count;++k) {
    auto line=refs[k];std::vector<double> le,ly,lm,ls,lb,lr,neighbor;double dev=0.;
    const double lo=k?std::max(guesses[k]-halfwidths[k],(guesses[k-1]+guesses[k])/2):edges.front();
    const double hi=k+1<count?std::min(guesses[k]+halfwidths[k],(guesses[k+1]+guesses[k])/2):edges.back();
    for(int i=0;i<n;++i)if((edges[i]+edges[i+1])/2>=lo&&(edges[i]+edges[i+1])/2<=hi) {
      if(le.empty())le.push_back(edges[i]);
      le.push_back(edges[i+1]);ly.push_back(counts[i]);lm.push_back(model[i]);ls.push_back(parts[k+1][i]);lb.push_back(parts[0][i]);neighbor.push_back(model[i]-parts[k+1][i]-parts[0][i]);
      const double m=std::max(model[i],1e-200),y=counts[i];const double d=std::max(0.,2*(m-y+(y>0?y*std::log(y/m):0.)));dev+=d;lr.push_back((y>=m?1:-1)*std::sqrt(d));
    }
    bool shape_bound=false;for(int z=0;z<stride;++z)if(bound.count(stride*k+z))shape_bound=true;
    const double area=std::exp(p[stride*k]),area_se=area*errors[stride*k],mu=p[stride*k+1],sigma=std::exp(p[stride*k+2]);
    const bool usable=best.at("success").get<bool>()&&best.at("covariance_status").get<int>()>=2&&!shape_bound&&
      errors[stride*k+1]>0&&errors[stride*k+1]<sigma&&area>5*area_se;
    const int free=stride-int(fixed.count(stride*k+3))-int(fixed.count(stride*k+4))-int(fixed.count(stride*k+5));const int dof=int(ly.size())-free-2;
    const double quality=dof>0?dev/dof:1e100;
    line.update({{"status","joint_peak_fit"},{"center",mu},{"center_se",errors[stride*k+1]},{"sigma",sigma},
      {"area",area},{"area_se",area_se},{"tail_model",dual?(p[stride*k+3]==0?"pure_left_tail":"double_left_tail"):(p[stride*k+3]==0?"gaussian":p[stride*k+3]==1?"pure_left_tail":"left_tail")},
      {"fraction",p[stride*k+3]},{"tau",std::exp(p[stride*k+4])},{"tau2",dual?json(std::exp(p[stride*k+5])):json(nullptr)},{"shape_bound",shape_bound},{"center_usable",usable},
      {"deviance",dev},{"ndof",dof},{"deviance_per_dof",quality},{"counts",ly},{"edges",le},{"model",lm},
      {"signal",ls},{"background",lb},{"neighbor_signal",neighbor},{"residual",lr},
      {"components_adc_offset_and_weight",components[k]},
      {"old_residual_keV",energy(mu,seed)-line.at("energy_keV").get<double>()},
      {"error_scope","conditional joint-fit covariance includes other source peaks; no model, source or dead-layer systematic uncertainty"},
      {"local_deviance_note","Descriptive local D/(bins - feature free parameters - 2); shared fit, not an independent chi-square test"}});
    const bool used=usable&&quality<=cfg.at("maximum_peak_deviance_per_dof").get<double>();line["used_for_calibration"]=used;
    if(used)anchors.push_back(line);else line["anchor_rejection_reason"]=shape_bound?"shape boundary":!usable?"fit/covariance/significance":"local residual structure";
    result["peaks"].push_back(line);
  }
  if(anchors.size()>=2) {
    auto curve=polynomial(anchors,1);const auto coeff=curve.at("coefficients").get<Coeff>();
    // Propagate the correlated source centers through the same equal-weight
    // line fit. A shared joint fit does not produce independent centers.
    const int m=int(anchors.size());TMatrixD design(m,2);TVectorD target(m);
    for(int i=0;i<m;++i){design(i,0)=1.;design(i,1)=anchors[i].at("center");target[i]=anchors[i].at("energy_keV");}
    TMatrixD transpose(TMatrixD::kTransposed,design),normal=transpose*design;Bool_t ok=false;TDecompSVD solver(normal);const auto inverse=solver.Invert(ok);require(ok,"Joint calibration covariance inversion failed");
    TMatrixD jacobian(2,m);jacobian.Zero();
    for(int i=0;i<m;++i){const double x=anchors[i].at("center"),r=anchors[i].at("energy_keV").get<double>()-energy(x,coeff);
      for(int a=0;a<2;++a)jacobian(a,i)=inverse(a,0)*(-coeff[1])+inverse(a,1)*(r-coeff[1]*x);}
    json covariance=json::array();
    for(int a=0;a<3;++a){json row=json::array();for(int b=0;b<3;++b){double value=0.;if(a<2&&b<2)for(int i=0;i<m;++i)for(int j=0;j<m;++j){int ki=-1,kj=-1;for(int k=0;k<count;++k){if(refs[k].at("id")==anchors[i].at("id"))ki=k;if(refs[k].at("id")==anchors[j].at("id"))kj=k;}value+=jacobian(a,i)*best.at("covariance")[stride*ki+1][stride*kj+1].get<double>()*jacobian(b,j);}row.push_back(value);}covariance.push_back(row);}
    curve["covariance"]=covariance;curve["weights"]="equal anchors; covariance propagates correlated joint-fit centers";
    curve["conditional_center_chi2_note"]="Legacy diagonal diagnostic only; source-center correlations are stored in joint_fit.covariance";
    const bool eligible=curve.at("max_residual_keV").get<double>()<=cfg.at("maximum_energy_residual_keV").get<double>()&&coeff[1]>0;
    result["calibration"]=curve;result["candidate_eligible"]=eligible;
    result["status"]=eligible?(anchors.size()==2?"two_point_candidate_unvalidated":"candidate"):"calibration_review_required";
    for(auto& line:result["peaks"])line["new_residual_keV"]=energy(line.at("center"),coeff)-line.at("energy_keV").get<double>();
  }
  return result;
}
}}
#endif
