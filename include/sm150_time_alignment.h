#ifndef SM150_TIME_ALIGNMENT_H
#define SM150_TIME_ALIGNMENT_H

// Versioned complete (one-pass) or historical additive timing parameters.
// Both CFD branches are explicit. No ROOT API is used; the JSON header is
// supplied by the same ROOT include path used by the production build.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <nlohmann/json.hpp>

namespace sm150 {
class TimeAlignment {
 public:
  TimeAlignment() : active_(false), total_(false), run_(0), s3_shift_(0.) {}

  void Load(const std::string& pars, int run) {
    std::ifstream input((pars + "/time_alignment/active.txt").c_str());
    if (!input) { Reset(run); return; }
    std::string profile, extra;
    if (!(input >> profile) || (input >> extra))
      throw std::runtime_error("Invalid time_alignment/active.txt");
    if (profile == "legacy") { Reset(run); return; }
    if (profile.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_") != std::string::npos)
      throw std::runtime_error("Invalid time alignment profile name");
    LoadProfile(pars + "/time_alignment/versions/" + profile, run);
    if (active_) profile_ = profile;
  }

  void LoadProfile(const std::string& directory, int run) {
    Reset(run);
    std::ifstream mode((directory + "/correction_mode.txt").c_str());
    if (mode) {
      std::string value, extra;
      if (!(mode >> value) || (mode >> extra) || value != "total")
        throw std::runtime_error("Invalid timing correction mode");
      total_ = true;
    }
    std::ifstream storage((directory + "/storage_format.txt").c_str());
    if (storage) {
      std::string value,extra;
      if (!(storage>>value) || (storage>>extra) || value!="json")
        throw std::runtime_error("Invalid timing storage format");
      LoadCompleteJson(directory,run);
      return;
    }
    std::ifstream references((directory + "/reference_offsets.tsv").c_str());
    if (!references) throw std::runtime_error("Missing reference_offsets.tsv: " + directory);
    std::string line, group;
    bool found = false;
    while (std::getline(references, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream row(line);
      int r; double shift;
      if (!(row >> r >> group >> shift) || !std::isfinite(shift))
        throw std::runtime_error("Invalid reference offset row");
      if (r == run) {
        if (found) throw std::runtime_error("Duplicate run in reference offsets");
        s3_shift_ = shift; group_ = group; found = true;
      }
    }
    // Runs outside the selected profile retain legacy processing. Never
    // interpolate parameters between run numbers or acquisition periods.
    if (!found) { Reset(run); return; }
    if (group_ != "high1" && group_ != "high2" && group_ != "low")
      throw std::runtime_error("Unsupported time alignment period");
    std::ostringstream run_name; run_name.width(4); run_name.fill('0'); run_name << run;
    std::ifstream data((directory + "/runs/" + run_name.str() + ".tsv").c_str());
    if (!data) throw std::runtime_error("Missing per-run time alignment table");
    while (std::getline(data, line)) {
      if (line.empty() || line[0] == '#') continue;
      int ge, state, enabled; double energy, delta;
      std::istringstream row(line);
      if (!(row >> ge >> state >> energy >> delta >> enabled) ||
          state < 0 || state > 1 || enabled < 0 || enabled > 1 ||
          !std::isfinite(energy) || !std::isfinite(delta) || energy < 40 || energy > 2000)
        throw std::runtime_error("Invalid Ge energy offset row");
      Curve& curve = curves_[std::make_pair(ge, state)];
      if (!curve.energy.empty() && (energy <= curve.energy.back() || curve.enabled != (enabled != 0)))
        throw std::runtime_error("Unordered or inconsistent Ge energy offset table");
      curve.enabled = enabled != 0;
      curve.energy.push_back(energy); curve.log_energy.push_back(std::log(energy)); curve.delta.push_back(delta);
    }
    if (curves_.size() != 40) throw std::runtime_error("Expected 20 Ge x 2 CFD branches");
    for (std::map<std::pair<int,int>,Curve>::const_iterator it=curves_.begin(); it!=curves_.end(); ++it) {
      if (it->second.energy.size()<2 || it->second.energy.front()!=40. || it->second.energy.back()!=2000.)
        throw std::runtime_error("Ge offset table lacks domain endpoints");
    }
    if (group_ == "low" || total_) {
      std::ifstream particles((directory + "/particles/" + run_name.str() + ".tsv").c_str());
      if (!particles) throw std::runtime_error("Missing low-beam particle-reference table");
      while (std::getline(particles,line)) {
        if (line.empty() || line[0]=='#') continue;
        std::istringstream row(line); int key; Piece piece;
        if (!(row>>key>>piece.lo>>piece.hi>>piece.slope>>piece.intercept) ||
            key<10000 || !std::isfinite(piece.lo) || !std::isfinite(piece.hi) ||
            !std::isfinite(piece.slope) || !std::isfinite(piece.intercept) ||
            piece.lo<0. || piece.hi<=piece.lo)
          throw std::runtime_error("Invalid low-beam particle-reference row");
        std::vector<Piece>& curve=particle_curves_[key];
        if ((!curve.empty() && piece.lo!=curve.back().hi) || (curve.empty() && piece.lo!=0.))
          throw std::runtime_error("Unordered low-beam particle-reference intervals");
        curve.push_back(piece);
      }
      if (particle_curves_.size()!=152) throw std::runtime_error("Expected all152 particle channels");
      for (std::map<int,std::vector<Piece> >::const_iterator it=particle_curves_.begin();it!=particle_curves_.end();++it)
        if (it->second.back().hi!=1e9) throw std::runtime_error("Particle-reference intervals lack upper bound");
      if (group_ == "low" || (total_ && ParticleEnergyDependent())) {
      std::ifstream energy((directory + "/particle_energy/" + run_name.str() + ".tsv").c_str());
      if (!energy) throw std::runtime_error("Missing frozen particle-energy calibration");
      while (std::getline(energy,line)) {
        if (line.empty() || line[0]=='#') continue;
        std::istringstream row(line); int key; std::vector<double> coefficients(3);
        if (!(row>>key>>coefficients[0]>>coefficients[1]>>coefficients[2]) ||
            particle_calibration_.count(key) || !particle_curves_.count(key) ||
            !std::isfinite(coefficients[0]) || !std::isfinite(coefficients[1]) || !std::isfinite(coefficients[2]))
          throw std::runtime_error("Invalid frozen particle-energy calibration row");
        particle_calibration_[key]=coefficients;
      }
      if (particle_calibration_.size()!=particle_curves_.size())
        throw std::runtime_error("Incomplete frozen particle-energy calibration");
      }
    }
    std::ifstream native((directory + (total_ ? "/channel_offsets_" : "/native_offsets_") + group_ + ".tsv").c_str());
    if (!native) throw std::runtime_error("Missing frozen native offset baseline");
    while (std::getline(native,line)) {
      if (line.empty() || line[0]=='#') continue;
      std::istringstream row(line); int key; double offset;
      if (!(row>>key>>offset) || !std::isfinite(offset) || native_.count(key))
        throw std::runtime_error("Invalid frozen native offset baseline");
      native_[key]=offset;
    }
    if (native_.empty()) throw std::runtime_error("Empty frozen native offset baseline");
    std::ifstream hash((directory + "/source_sha256.txt").c_str());
    if (!(hash>>source_sha_) || source_sha_.size()!=64)
      throw std::runtime_error("Missing calibration source identity");
    profile_ = directory.substr(directory.find_last_of('/')+1);
    active_ = true;
  }

  void ValidateNativeOffsets(const std::map<int,double>& actual) const {
    if (!active_ || total_) return;
    for (std::map<int,double>::const_iterator it=native_.begin(); it!=native_.end(); ++it) {
      std::map<int,double>::const_iterator found=actual.find(it->first);
      if (found==actual.end() || !std::isfinite(found->second) || std::abs(found->second-it->second)>1e-8)
        throw std::runtime_error("Native offsets changed since calibration; do not add the new profile twice");
    }
  }

  bool GeDelta(int key, bool failed, double energy_keV, double& delta_ns) const {
    if (total_) throw std::runtime_error("Complete profile requires TotalGeOffset, not an added correction");
    delta_ns=0.;
    if (!active_ || !std::isfinite(energy_keV) || energy_keV<40. || energy_keV>=2000.) return false;
    std::map<std::pair<int,int>,Curve>::const_iterator it=curves_.find(std::make_pair(key,failed?1:0));
    if (it==curves_.end() || !it->second.enabled) return false;
    const Curve& curve=it->second;
    std::vector<double>::const_iterator hi=std::upper_bound(curve.energy.begin(),curve.energy.end(),energy_keV);
    const std::size_t j=static_cast<std::size_t>(hi-curve.energy.begin()-1);
    const double fraction=(std::log(energy_keV)-curve.log_energy[j])/(curve.log_energy[j+1]-curve.log_energy[j]);
    delta_ns=curve.delta[j]+fraction*(curve.delta[j+1]-curve.delta[j]);
    return true;
  }

  // One common shift for all S3 rings AND sectors. Never shift just ring1-8:
  // that would introduce an artificial front/back time difference.
  double ParticleDelta(bool is_s3) const {
    if (total_) throw std::runtime_error("Complete profile requires a complete particle offset");
    return active_ && is_s3 ? -s3_shift_ : 0.;
  }
  // Low-beam candidate references depend on physical particle energy and channel.
  // These are inherited, separately recorded terms, not later pair-level fit
  // nuisances. High-beam callers keep exactly the existing common S3 shift.
  double ParticleDelta(int key, double energy_keV, bool is_s3) const {
    if (total_) throw std::runtime_error("Complete profile requires TotalParticleOffset");
    if (!active_ || group_!="low") return ParticleDelta(is_s3);
    if (!std::isfinite(energy_keV)) throw std::runtime_error("Non-finite particle timing energy");
    std::map<int,std::vector<Piece> >::const_iterator found=particle_curves_.find(key);
    if (found==particle_curves_.end()) throw std::runtime_error("Unmapped low-beam particle channel");
    const std::vector<Piece>& pieces=found->second;
    if (energy_keV>=1e9) throw std::runtime_error("Particle timing energy outside calibrated ADC domain");
    // Nonpositive baseline/noise energies are outside particle-energy curves;
    // retain their channel/run constant, as in the frozen source evaluator.
    energy_keV=std::max(0.,energy_keV);
    for (std::size_t i=0;i<pieces.size();++i)
      if (energy_keV>=pieces[i].lo && energy_keV<pieces[i].hi)
        return pieces[i].slope*energy_keV+pieces[i].intercept;
    throw std::runtime_error("Missing low-beam particle timing interval");
  }
  double ParticleDeltaFromAdc(int key, double adc, bool is_s3) const {
    if (total_) throw std::runtime_error("Complete profile requires TotalParticleOffsetFromAdc");
    if (!active_ || group_!="low") return ParticleDelta(is_s3);
    std::map<int,std::vector<double> >::const_iterator found=particle_calibration_.find(key);
    if (found==particle_calibration_.end()) throw std::runtime_error("Unmapped low-beam particle calibration");
    const std::vector<double>& c=found->second;
    return ParticleDelta(key,c[0]+c[1]*adc+c[2]*adc*adc,is_s3);
  }
  // These return already merged, complete parameters. Never add an old offset.
  double TotalGeOffset(int key, bool failed, double energy_keV) const {
    RequireTotal();
    if (!std::isfinite(energy_keV) || energy_keV<40. || energy_keV>=2000.) return ChannelFallback(key);
    std::map<std::pair<int,int>,Curve>::const_iterator it=curves_.find(std::make_pair(key,failed?1:0));
    if (it==curves_.end()) return ChannelFallback(key);
    const Curve& curve=it->second;
    // Held branches also contain their full frozen constant in this table.
    const std::size_t j=static_cast<std::size_t>(std::upper_bound(curve.energy.begin(),curve.energy.end(),energy_keV)-curve.energy.begin()-1);
    const double f=(std::log(energy_keV)-curve.log_energy[j])/(curve.log_energy[j+1]-curve.log_energy[j]);
    double value=curve.delta[j]+f*(curve.delta[j+1]-curve.delta[j]);
    if (failed) {
      std::map<int,CenterModel>::const_iterator model=center_models_.find(key);
      if (model!=center_models_.end()) value+=EvaluateCenter(model->second,energy_keV);
    }
    return value;
  }
  double TotalParticleOffset(int key, double energy_keV) const {
    RequireTotal();
    if (!std::isfinite(energy_keV) || energy_keV>=1e9)
      throw std::runtime_error("Particle timing energy outside calibrated ADC domain");
    std::map<int,std::vector<Piece> >::const_iterator it=particle_curves_.find(key);
    if (it==particle_curves_.end()) throw std::runtime_error("Unmapped particle timing channel");
    energy_keV=std::max(0.,energy_keV);
    for (std::size_t i=0;i<it->second.size();++i) {
      const Piece& p=it->second[i];
      if (energy_keV>=p.lo && energy_keV<p.hi) return p.slope*energy_keV+p.intercept;
    }
    throw std::runtime_error("Missing complete particle timing interval");
  }
  double TotalParticleOffsetFromAdc(int key, double adc) const {
    RequireTotal();
    if (!std::isfinite(adc)) throw std::runtime_error("Non-finite particle ADC");
    std::map<int,std::vector<double> >::const_iterator it=particle_calibration_.find(key);
    if (it==particle_calibration_.end()) {
      // Support legacy constant high-beam profiles without energy-calibration tables; energy-dependent profiles must not fall back to E=0.
      std::map<int,std::vector<Piece> >::const_iterator p=particle_curves_.find(key);
      if (group_=="low" || p==particle_curves_.end() || p->second.empty())
        throw std::runtime_error("Unmapped physical particle calibration");
      for (const Piece& piece:p->second)
        if (piece.slope!=0. || piece.intercept!=p->second.front().intercept)
          throw std::runtime_error("Energy-dependent particle timing requires ADC calibration");
      return TotalParticleOffset(key,0.);
    }
    const std::vector<double>& c=it->second;
    return TotalParticleOffset(key,c[0]+c[1]*adc+c[2]*adc*adc);
  }
  bool total() const { return active_ && total_; }
  bool active() const { return active_; }
  const std::string& profile() const { return profile_; }
  const std::string& source_sha256() const { return source_sha_; }

 private:
  struct Curve { bool enabled=false; std::vector<double> energy,log_energy,delta; };
  struct Piece { double lo,hi,slope,intercept; };
  struct CenterModel { std::string family; double lo,hi; std::vector<double> coefficients,knots; };
  static double EvaluateCenter(const CenterModel& m,double energy) {
    const double e=std::max(m.lo,std::min(m.hi,energy));
    const double u=std::log(e/334.); const std::vector<double>& c=m.coefficients;
    if (m.family=="constant") return c[0];
    if (m.family=="log_linear") return c[0]+c[1]*u;
    if (m.family=="log_quadratic") return c[0]+c[1]*u+c[2]*u*u;
    if (m.family=="inverse_sqrt") return c[0]+c[1]*(std::sqrt(334./e)-1.);
    if (m.family=="inverse") return c[0]+c[1]*(334./e-1.);
    if (m.family=="piecewise_log_linear") {
      if (e<=m.knots.front()) return c.front();
      if (e>=m.knots.back()) return c.back();
      const std::size_t j=std::upper_bound(m.knots.begin(),m.knots.end(),e)-m.knots.begin()-1;
      const double f=std::log(e/m.knots[j])/std::log(m.knots[j+1]/m.knots[j]);
      return c[j]+f*(c[j+1]-c[j]);
    }
    throw std::runtime_error("Unsupported failed-CFD center model");
  }
  void LoadCompleteJson(const std::string& directory,int run) {
    if (!total_) throw std::runtime_error("JSON timing profile must be complete");
    std::ifstream input((directory+"/calibration.json").c_str());
    if (!input) throw std::runtime_error("Missing complete timing JSON");
    nlohmann::json data; input>>data;
    if (data.at("correction_mode")!="total" || data.at("time_definition")!="t_particle - t_Ge")
      throw std::runtime_error("Unsupported JSON timing convention");
    std::ostringstream name;name.width(4);name.fill('0');name<<run;
    const std::string r=name.str();
    if (!data.at("groups").contains(r)) { Reset(run);return; }
    group_=data.at("groups").at(r).get<std::string>();
    if (group_!="low" && group_!="high1" && group_!="high2") throw std::runtime_error("Invalid timing period");
    const std::vector<int> ge=data.at("ge").get<std::vector<int> >();
    if (ge.size()!=20) throw std::runtime_error("Expected20 Ge channels");
    for (std::size_t i=0;i<ge.size();++i) for (int state=0;state<2;++state) {
      const std::string key=r+":"+std::to_string(ge[i])+":"+std::to_string(state);
      const nlohmann::json& source=data.at("ge_curves").at(key);Curve curve;
      curve.energy=source.at("knots_keV").get<std::vector<double> >();
      curve.delta=source.at("offset_ns").get<std::vector<double> >();
      curve.enabled=true;
      if (curve.energy.size()<2 || curve.delta.size()!=curve.energy.size() || curve.energy.front()!=40. || curve.energy.back()!=2000.)
        throw std::runtime_error("Invalid complete Ge JSON curve");
      for (std::size_t j=0;j<curve.energy.size();++j) {
        if (!std::isfinite(curve.delta[j]) || !std::isfinite(curve.energy[j]) || (j && curve.energy[j]<=curve.energy[j-1]))
          throw std::runtime_error("Unordered/nonfinite Ge JSON curve");
        curve.log_energy.push_back(std::log(curve.energy[j]));
      }
      curves_[std::make_pair(ge[i],state)]=curve;
    }
    const nlohmann::json& particles=data.at("particle_tables").at(r);
    for (nlohmann::json::const_iterator it=particles.begin();it!=particles.end();++it) {
      const int key=std::stoi(it.key());std::vector<Piece>& pieces=particle_curves_[key];
      for (std::size_t j=0;j<it.value().size();++j) {
        const std::vector<double> v=it.value().at(j).get<std::vector<double> >();
        if (v.size()!=4 || !std::all_of(v.begin(),v.end(),[](double x){return std::isfinite(x);}) || v[0]<0 || v[1]<=v[0] || (j==0 ? v[0]!=0. : v[0]!=pieces.back().hi))
          throw std::runtime_error("Invalid JSON particle intervals");
        Piece p={v[0],v[1],v[2],v[3]};pieces.push_back(p);
      }
      if (pieces.empty() || pieces.back().hi!=1e9) throw std::runtime_error("Incomplete JSON particle domain");
    }
    if (particle_curves_.size()!=152) throw std::runtime_error("Expected152 particle channels");
    // Allow complete energy corrections for both beam intensities while retaining support for legacy constant high-beam JSON profiles.
    const bool has_particle_calibration=data.contains("particle_calibration") &&
      data.at("particle_calibration").contains(r);
    if (!has_particle_calibration && (group_=="low" || ParticleEnergyDependent()))
      throw std::runtime_error("Missing JSON calibration for energy-dependent particle timing");
    if (has_particle_calibration) {
      const nlohmann::json& cal=data.at("particle_calibration").at(r);
      for (nlohmann::json::const_iterator it=cal.begin();it!=cal.end();++it) {
        const int key=std::stoi(it.key());const std::vector<double> c=it.value().get<std::vector<double> >();
        if (c.size()!=3 || !particle_curves_.count(key) || !std::all_of(c.begin(),c.end(),[](double x){return std::isfinite(x);}))
          throw std::runtime_error("Invalid JSON physical particle calibration");
        particle_calibration_[key]=c;
      }
      if (particle_calibration_.size()!=152) throw std::runtime_error("Incomplete JSON particle calibration");
    }
    const nlohmann::json& constants=data.at("channel_offsets").at(group_);
    for (nlohmann::json::const_iterator it=constants.begin();it!=constants.end();++it) {
      const double v=it.value().get<double>();if (!std::isfinite(v)) throw std::runtime_error("Nonfinite channel offset");
      native_[std::stoi(it.key())]=v;
    }
    if (data.contains("failed_center_models")) for (std::size_t i=0;i<ge.size();++i) {
      const std::string run_key=r+":"+std::to_string(ge[i]);
      const std::string period_key=group_+":"+std::to_string(ge[i]);
      const std::string key=data.at("failed_center_models").contains(run_key) ? run_key : period_key;
      if (!data.at("failed_center_models").contains(key)) continue;
      const nlohmann::json& source=data.at("failed_center_models").at(key);CenterModel m;
      m.family=source.at("family").get<std::string>();m.lo=source.at("support_keV").at(0).get<double>();m.hi=source.at("support_keV").at(1).get<double>();
      m.coefficients=source.at("coefficients_ns").get<std::vector<double> >();
      const std::size_t n=m.coefficients.size();
      if (!(40.<=m.lo && m.lo<m.hi && m.hi<=2000.) || !std::all_of(m.coefficients.begin(),m.coefficients.end(),[](double x){return std::isfinite(x);}))
        throw std::runtime_error("Invalid failed-center coefficients/support");
      if (m.family=="piecewise_log_linear") {
        m.knots=source.at("knots_keV").get<std::vector<double> >();
        if (n<2 || m.knots.size()!=n || m.knots.front()!=m.lo || m.knots.back()!=m.hi)
          throw std::runtime_error("Invalid failed-center knots");
        for (std::size_t j=0;j<n;++j) if (!std::isfinite(m.knots[j]) || (j && m.knots[j]<=m.knots[j-1]))
          throw std::runtime_error("Unordered failed-center knots");
      } else if (!((m.family=="constant" && n==1) || ((m.family=="log_linear" || m.family=="inverse" || m.family=="inverse_sqrt") && n==2) || (m.family=="log_quadratic" && n==3)))
        throw std::runtime_error("Unsupported failed-center function/count");
      center_models_[ge[i]]=m;
    }
    std::ifstream hash((directory+"/source_sha256.txt").c_str());
    if (!(hash>>source_sha_) || source_sha_.size()!=64) throw std::runtime_error("Missing JSON calibration identity");
    profile_=directory.substr(directory.find_last_of('/')+1);active_=true;
  }
  void RequireTotal() const {
    if (!active_ || !total_) throw std::runtime_error("A complete active timing profile is required");
  }
  bool ParticleEnergyDependent() const {
    for (const auto& channel:particle_curves_)
      for (const Piece& piece:channel.second)
        if (piece.slope!=0. || piece.intercept!=channel.second.front().intercept) return true;
    return false;
  }
  double ChannelFallback(int key) const {
    std::map<int,double>::const_iterator it=native_.find(key);
    if (it==native_.end()) throw std::runtime_error("Unmapped complete channel timing constant");
    return it->second;
  }
  void Reset(int run) {
    active_=false;total_=false;run_=run;s3_shift_=0.;curves_.clear();native_.clear();
    particle_curves_.clear();particle_calibration_.clear();center_models_.clear();
    profile_="legacy";source_sha_.clear();group_.clear();
  }
  bool active_,total_; int run_; double s3_shift_;
  std::string profile_,source_sha_,group_;
  std::map<std::pair<int,int>,Curve> curves_;
  std::map<int,double> native_;
  std::map<int,std::vector<Piece> > particle_curves_;
  std::map<int,std::vector<double> > particle_calibration_;
  std::map<int,CenterModel> center_models_;
};
} // namespace sm150
#endif
