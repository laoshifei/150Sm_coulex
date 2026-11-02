// Joint residual gain for the two faces of a matched S3 pixel.
// Included after the common signed-spectrum implementation in conditioned_calibrate.cpp.
#ifndef SM150_S3_PIXEL_FIT_H
#define SM150_S3_PIXEL_FIT_H
namespace {
struct PixelComparison {
  std::array<Comparison,2> faces;
  double operator()(double q)const {
    const double a=faces[0](q),b=faces[1](q);
    return std::sqrt((a*a+b*b)/2.);
  }
  bool usable()const{return faces[0].usable()&&faces[1].usable();}
};
Estimate pixel_estimate(const PixelComparison& p) {
  return estimate([&](double q){return p(q);},.5,3.,181);
}
json pixel_metrics(const PixelComparison& p,double base,double candidate) {
  json f=json::array();bool supported=true;
  for(const auto& c:p.faces) {
    const auto r=c.counts(true),b=c.counts(false,base),n=c.counts(false,candidate);
    supported&=r[0]>=15&&b[0]>=15&&n[0]>=15;
    f.push_back({{"baseline_D",c(base)},{"candidate_D",c(candidate)},
      {"reference_counts",r},{"baseline_counts",b},{"candidate_counts",n}});
  }
  return {{"baseline_D",p(base)},{"candidate_D",p(candidate)},{"faces",f},{"supported",supported}};
}
std::vector<std::string> pixel_guard(const json& m,bool selection) {
  std::vector<std::string> flags;
  if(!m.at("supported").get<bool>())flags.push_back("insufficient_comparison_counts");
  const double b=m.at("baseline_D"),n=m.at("candidate_D");
  if(selection&&(n>.95*b||b-n<.001))flags.push_back("no_clear_selection_improvement");
  if(!selection&&n>b+.002)flags.push_back("confirmation_worse");
  for(size_t face=0;face<2;++face)
    if(double(m.at("faces")[face].at("candidate_D"))>double(m.at("faces")[face].at("baseline_D"))+.002)
      flags.push_back(face?"sector_worse":"ring_worse");
  return flags;
}
void pixel_self_test() {
  for(double truth:{.72,1.,1.618,2.35}) {
    PixelComparison p;
    for(int face=0;face<2;++face) {
      auto& c=p.faces[face];c.fixed_normalization=true;c.ar=face?250.:-150.;c.br=face?1.12:.94;
      c.ac=face?-320.:170.;c.bc=face?.91:1.07;
      for(int i=0;i<=bins;++i) {
        const double e=i*width;
        auto integral=[](double x){return 50000.*.5*(1.+std::erf((x-60000.)/(2500.*std::sqrt(2.))));};
        c.reference.r[i]=double(i)*2.;c.current.r[i]=double(i)*3.;
        c.reference.p[i]=integral(c.ar+c.br*e)+alpha*c.reference.r[i];
        c.current.p[i]=integral(truth*(c.ac+c.bc*e))+alpha*c.current.r[i];
      }
    }
    auto fit=pixel_estimate(p);
    require(std::abs(fit.q/truth-1.)<.001,"Joint pixel gain recovery failed");
    require(p.usable(),"Synthetic pixel unexpectedly unsupported");
    std::cout<<"pixel-self-test truth="<<truth<<" recovered="<<fit.q<<" distance="<<fit.value<<'\n';
  }
  std::cout<<"Joint pixel synthetic checks passed"<<std::endl;
}
void fit_pixels(const Data& current,const Data& ref,const json& package,
                const fs::path& input,const fs::path& reference,const fs::path& baseline,
                const fs::path& output,int first_ring,int last_ring) {
  const auto start=std::chrono::steady_clock::now();
  const std::string run=current.meta.at("run"),ref_run=ref.meta.at("run");
  require(current.meta.at("period")=="low"&&ref.meta.at("period")=="low","Pixel study expects low beam");
  const auto& m=package.at("runs").at(run);const auto& rm=package.at("runs").at(ref_run);
  require(m.at("active_rings")==24&&rm.at("active_rings")==24,"Incomplete low-beam ring coverage");
  json rows=json::array();int accepted=0;
  for(int ring=first_ring-1;ring<last_ring;++ring)for(int sector=0;sector<32;++sector) {
    const int pixel=ring*32+sector;
    auto make=[&](int fold) {
      PixelComparison p;
      for(int face=0;face<2;++face) {
        const int index=face?24+sector:ring;
        auto& c=p.faces[face];c.fixed_normalization=true;c.current=combine(current,face,pixel,fold);c.reference=combine(ref,face,pixel,fold);
        c.ac=m.at("A_keV")[index];c.bc=m.at("B")[index];
        const double rq=rm.at("q_position")[ring];c.ar=rq*double(rm.at("A_keV")[index]);c.br=rq*double(rm.at("B")[index]);
      }
      return p;
    };
    const auto train=make(-1);const double base=m.at("q_position")[ring];
    json row={{"run",run},{"reference_run",ref_run},{"ring",ring+1},{"sector",sector+1},
      {"pixel",pixel},{"baseline_q",base},{"q_selected",base},{"training_usable",train.usable()}};
    if(run==ref_run) {
      row.update({{"q_fitted",base},{"status","fixed_reference"},{"flags",json::array()}});rows.push_back(row);continue;
    }
    if(!train.usable()) {
      row.update({{"q_fitted",nullptr},{"status","unresolved_inherited"},{"flags",{"insufficient_training_counts"}}});
      rows.push_back(row);continue;
    }
    const auto fit=pixel_estimate(train),f0=pixel_estimate(make(0)),f2=pixel_estimate(make(2));
    auto alternate=train;for(auto& c:alternate.faces){c.lo=30000;c.hi=80000;}
    const auto alt=pixel_estimate(alternate);
    std::array<double,2> face_q{};
    for(int face=0;face<2;++face)face_q[face]=estimate([&](double q){return train.faces[face](q);},.5,3.,181).q;
    const double split=std::abs(std::log(f0.q/f2.q)),window=std::abs(std::log(alt.q/fit.q));
    const double face_delta=std::abs(std::log(face_q[0]/face_q[1])),support=std::log(fit.high/fit.low);
    std::vector<std::string> flags;
    if(fit.boundary||fit.q<.505||fit.q>2.97)flags.push_back("search_boundary");
    if(fit.low<=.500001||fit.high>=2.999999)flags.push_back("profile_reaches_search_boundary");
    if(split>.05)flags.push_back("fit_segment_difference_above5pct");
    if(window>.05)flags.push_back("energy_window_difference_above5pct");
    if(face_delta>.05)flags.push_back("face_gain_difference_above5pct");
    if(support>.20)flags.push_back("broad_descriptive_profile");
    const auto sel=pixel_metrics(make(1),base,fit.q),check=pixel_metrics(make(3),base,fit.q);
    const auto sf=pixel_guard(sel,true),cf=pixel_guard(check,false);
    flags.insert(flags.end(),sf.begin(),sf.end());flags.insert(flags.end(),cf.begin(),cf.end());
    const bool use=flags.empty();if(use)++accepted;
    row.update({{"q_fitted",fit.q},{"q_selected",use?fit.q:base},{"status",use?"pixel_refit_accepted":"unresolved_inherited"},
      {"flags",flags},{"q_fit_fold0",f0.q},{"q_fit_fold2",f2.q},{"q_alternate_window",alt.q},{"q_each_face",face_q},
      {"split_log_difference",split},{"window_log_difference",window},{"face_log_difference",face_delta},
      {"profile_support",{fit.low,fit.high}},{"profile",{{"q",fit.grid},{"D",fit.profile}}},
      {"training",pixel_metrics(train,base,fit.q)},{"selection",sel},{"confirmation",check}});
    rows.push_back(row);
  }
  json result={{"schema","sm150_s3_joint_pixel_fit_v1"},{"status","complete"},{"run",run},{"reference_run",ref_run},
    {"rings",{first_ring,last_ring}},{"pixels",rows},{"accepted_pixels",accepted},
    {"input",fs::absolute(input).string()},{"reference",fs::absolute(reference).string()},
    {"input_sha256",hash(read_text(input))},{"reference_sha256",hash(read_text(reference))},{"baseline_sha256",hash(read_text(baseline))},
    {"formula","Er=q[ring,sector]*(Ar+Br*E0r); Es=q[ring,sector]*(As+Bs*E0s)"},
    {"objective","Equal face mean of squared signed cumulative-spectrum distances on20-90MeV; fixed full E0 histogram net normalization, P-0.5R"},
    {"q_bounds",{.5,3.}},{"validation_status","Historical folds; fresh disjoint segments must be checked separately"},
    {"program_sha256",hash(read_text(fs::canonical("/proc/self/exe")))},
    {"elapsed_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()}};
  fs::create_directories(fs::absolute(output).parent_path());const fs::path partial=output.string()+".partial";
  std::ofstream out(partial);out<<result.dump(2)<<'\n';out.close();require(bool(out),"Pixel fit write failed");
  require(::link(partial.c_str(),output.c_str())==0&&::unlink(partial.c_str())==0,"Pixel output publication failed");
  std::cout<<run<<" joint pixel fit: "<<accepted<<"/"<<rows.size()<<" updated; seconds="<<result.at("elapsed_seconds")<<std::endl;
}
}
#endif
