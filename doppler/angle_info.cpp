#include "angle_info.h"
#include "../include/sm150_parameter_paths.h"
#include <iostream>
#include <fstream>

#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "set.h"

angle_info::angle_info()
{
  if(!InitAngleInfo()){
    throw std::invalid_argument("can not init angle info.");
  }
  PrintAngleInfo();
  
  if(!InitEinEoutInfo()){
    throw std::invalid_argument("can not init e_in and e_out info.");
  }

  PrintEinEoutInfo();
  PrintThetaP2RInfo();

  PrintBetaInfo();
}

//
angle_info::~angle_info()
{
  
}

//  
bool angle_info::InitAngleInfo()
{
  std::cout << "start read angle info data" << std::endl;

  std::string line;
  int ring, sector;
  double theta, phi;
  int key = 0;

  //
  std::ifstream fi_ge;
  fi_ge.open("../pars/angle_info/ge_angle.txt");
  if(!fi_ge){
    std::cout << "can not open ge_angle.txt." << std::endl;
    return 0;
  }

  std::getline(fi_ge, line);
  while(1){
    fi_ge >> ring  >> sector >> theta >> phi;
    if(!fi_ge.good()) break;

    key = 100*ring+sector;
    map_ge_angle[key] = std::make_pair(theta, phi);
  }
  fi_ge.close();

  //
  std::ifstream fi_spider;
  fi_spider.open("../pars/angle_info/spider_angle.txt");
  if(!fi_spider){
    std::cout << "can not open spider_angle.txt." << std::endl;
    return 0;
  }

  std::getline(fi_spider, line);
  while(1){
    fi_spider >> ring  >> sector >> theta >> phi;
    if(!fi_spider.good()) break;

    key = 100*ring+sector;
    map_spider_angle[key] = std::make_pair(theta, phi);
  }
  fi_spider.close();

  //
  std::ifstream fi_s3_sector;
  fi_s3_sector.open("../pars/angle_info/s3_sector_angle.txt");
  if(!fi_s3_sector){
    std::cout << "can not open s3_sector_angle.txt." << std::endl;
    return 0;
  }

  std::getline(fi_s3_sector, line);
  while(1){
    fi_s3_sector >> sector >> phi;
    if(!fi_s3_sector.good()) break;

    key = sector;
    map_s3_sector_phi[key] = phi;
  }
  fi_s3_sector.close();

  //
  std::ifstream fi_s3_ring;
  fi_s3_ring.open("../pars/angle_info/s3_ring_angle.txt");
  if(!fi_s3_ring){
    std::cout << "can not open s3_ring_angle.txt." << std::endl;
    return 0;
  }

  std::getline(fi_s3_ring, line);
  while(1){
    fi_s3_ring >> ring >> theta;
    if(!fi_s3_ring.good()) break;

    key = ring;
    map_s3_ring_theta[key] = theta;
  }
  fi_s3_ring.close();

  return 1;
}

//
void angle_info::PrintAngleInfo()
{
  std::cout << "print angle info ..." << std::endl;
  
  //
  std::cout << std::endl;
  std::cout << "print ge angle info ..." << std::endl;
  std::map<int, std::pair<double, double>>::iterator it_ge = map_ge_angle.begin();
  for(;it_ge!=map_ge_angle.end();it_ge++){
    std::cout << "ring " << it_ge->first/100 << " sector " << it_ge->first%100<< " => " << it_ge->second.first << " " << it_ge->second.second << '\n';
  }

  //
  std::cout << std::endl;
  std::cout << "print spider angle info ..." << std::endl;
  std::map<int, std::pair<double, double>>::iterator it_spider = map_spider_angle.begin();
  for(;it_spider!=map_spider_angle.end();it_spider++){
    std::cout << "ring " << it_spider->first/100 << " sector " << it_spider->first%100<< " => " << it_spider->second.first << " " << it_spider->second.second << '\n';
  }

  //
  std::cout << std::endl;
  std::cout << "print s3 sector phi info ..." << std::endl;
  std::map<int, double>::iterator it_s3_sector = map_s3_sector_phi.begin();
  for(;it_s3_sector!=map_s3_sector_phi.end();it_s3_sector++){
    std::cout << "sector " << it_s3_sector->first << " => " << it_s3_sector->second << '\n';
  }

  //
  std::cout << std::endl;
  std::cout << "print s3 ring theta info ..." << std::endl;
  std::map<int, double>::iterator it_s3_ring = map_s3_ring_theta.begin();
  for(;it_s3_ring!=map_s3_ring_theta.end();it_s3_ring++){
    std::cout << "ring " << it_s3_ring->first << " => " << it_s3_ring->second << '\n';
  }
}

//
bool angle_info::InitEinEoutInfo()
{
  std::cout << "start read Ein Eout info data" << std::endl;

  // The positional interface keeps its frozen historical configuration. The
  // current build interface below reads precomputed beta values instead.
  std::ifstream legacy_file("../pars/doppler/angle_table_baseline.json");
  if(!legacy_file)return false;
  nlohmann::json legacy_parameters;legacy_file>>legacy_parameters;
  const double legacy_s3_scale=legacy_parameters.at("S3").at("energy_scale");
  if(!std::isfinite(legacy_s3_scale)||legacy_s3_scale<=0)return false;

  //
  std::ifstream fi_p;
  fi_p.open("../pars/beta_info/theta_ein_eout_p.txt");
  if(!fi_p){
    std::cout << "can not open theta_ein_eout_p.txt." << std::endl;
    return 0;
  }

  char c_theta_p[32];
  double theta_p, e_in, e_out;
  while(1){
    fi_p >> theta_p >> e_in >> e_out;
    if(!fi_p.good()) break;

    theta_p = std::round(theta_p * 100.0) / 100.0;
    std::sprintf(c_theta_p, "%.2f", theta_p);
    map_theta_ein_eout_p[c_theta_p] = std::make_pair(e_in, e_out);
  }
  fi_p.close();

  //
  std::ifstream fi_r;
  fi_r.open("../pars/beta_info/theta_ein_eout_r.txt");
  if(!fi_r){
    std::cout << "can not open theta_ein_eout_r.txt." << std::endl;
    return 0;
  }

  char c_theta_r[32];
  double theta_r;
  while(1){
    fi_r >> theta_p >> theta_r >> e_in >> e_out;
    if(!fi_r.good()) break;

    theta_p = std::round(theta_p * 100.0) / 100.0;
    theta_r = std::round(theta_r * 100.0) / 100.0;
    std::sprintf(c_theta_p, "%.2f", theta_p);
    std::sprintf(c_theta_r, "%.2f", theta_r);
    if(theta_p<90){
      map_theta_ein_eout_r[c_theta_r] = std::make_pair(e_in, e_out*legacy_s3_scale);
    }else{
      map_theta_ein_eout_r[c_theta_r] = std::make_pair(e_in, e_out);
    }
    map_theta_p2r[c_theta_p] = theta_r;
  }
  fi_r.close();

  return 1;
}

//
void angle_info::PrintEinEoutInfo()
{
  std::cout << "print ein eout info ..." << std::endl;
  
  //
  std::cout << std::endl;
  std::map<std::string, std::pair<double, double>>::iterator it_p = map_theta_ein_eout_p.begin();
  for(;it_p!=map_theta_ein_eout_p.end();it_p++){
    std::cout << it_p->first << " => " << it_p->second.first << " " << it_p->second.second << '\n';
  }

  //
  std::cout << std::endl;
  std::map<std::string, std::pair<double, double>>::iterator it_r = map_theta_ein_eout_r.begin();
  for(;it_r!=map_theta_ein_eout_r.end();it_r++){
    std::cout << it_r->first << " => " << it_r->second.first << " " << it_r->second.second << '\n';
  }
}

//
void angle_info::PrintThetaP2RInfo()
{
  std::cout << "print theta_p to theta_r info ..." << std::endl;

  //
  std::cout << std::endl;
  std::map<std::string, double>::iterator it = map_theta_p2r.begin();
  for(;it!=map_theta_p2r.end();it++){
    std::cout << it->first << " => " << it->second << '\n';
  }
}

//
void angle_info::PrintBetaInfo()
{
  std::cout << "print projectile beta info ..." << std::endl;

  //
  std::map<std::string, double>::iterator it = map_theta_p2r.begin();
  for(;it!=map_theta_p2r.end();it++){
    std::cout << "p ==> " << it->first << "deg, beta " << ComputeProjectileBeta(std::stod(it->first), 1000.) << std::endl;
    std::cout << "r ==> " << it->second << "deg, beta " << ComputeRecoilBeta(it->second, 1000.) << std::endl;
    std::cout << std::endl;
  }

  std::cout << std::endl;
}

//
double angle_info::GetRecoilTheta(double theta_p)
{
  char c_theta_p[32];
  theta_p = std::round(theta_p * 100.0) / 100.0;
  std::sprintf(c_theta_p, "%.2f", theta_p);

  return map_theta_p2r[c_theta_p];
}

//
double angle_info::GetRecoilPhi(double phi_p)
{ 
  return phi_p+180; 
}

// dt unit fs
double angle_info::ComputeProjectileBeta(double theta, double dt)
{
  double half_target_thickness = TARGETTHICKNESS/2.;

  char c_theta[32];
  theta = std::round(theta * 100.0) / 100.0;
  std::sprintf(c_theta, "%.2f", theta);

  double e_in = map_theta_ein_eout_p[c_theta].first;
  double e_out = map_theta_ein_eout_p[c_theta].second;

  double d = 0;
  if(theta>90)  d = half_target_thickness/std::cos((180.-theta)/180.*PI); // half target thickness = 0.8 um.
  else d = half_target_thickness/std::cos(theta/180.*PI); // half target thickness = 0.8 um.

  double m_p = 16*938.3 + 16*931.5 - 32*8.493;
  double beta_in = std::sqrt(2.*e_in/m_p);
  double beta_out = std::sqrt(2.*e_out/m_p);

  double v_in = beta_in*300; // d/v = time with unit ps
  double v_out = beta_out*300;

  double a = (v_out*v_out-v_in*v_in)/2./d;
  double v_now = v_in + a*dt/1000.;
  double beta_now = 0;
  
  if(v_now<v_out){
    beta_now = beta_out;
  }else{
    beta_now = v_now/300.;
  }

  // std::cout << "projectile a " << a << " t " << (v_out-v_in)/a << "ns " << v_now << " " << beta_now << std::endl;
  return beta_now;
}

// dt unit fs
double angle_info::ComputeRecoilBeta(double theta, double dt)
{
  double half_target_thickness = TARGETTHICKNESS/2.;

  char c_theta[32];
  theta = std::round(theta * 100.0) / 100.0;
  std::sprintf(c_theta, "%.2f", theta);

  double e_in = map_theta_ein_eout_r[c_theta].first;
  double e_out = map_theta_ein_eout_r[c_theta].second;

  double d = 0;
  if(theta>90)  d = half_target_thickness/std::cos((180.-theta)/180.*PI); // half target thickness = 0.8 um.
  else d = half_target_thickness/std::cos(theta/180.*PI); // half target thickness = 0.8 um.

  double m_r = 62*938.3 + 88*931.5 - 150*8.262;
  double beta_in = std::sqrt(2.*e_in/m_r);
  double beta_out = std::sqrt(2.*e_out/m_r);

  double v_in = beta_in*300; // d/v = time with unit ps
  double v_out = beta_out*300;

  double a = (v_out*v_out-v_in*v_in)/2./d;
  double v_now = v_in + a*dt/1000.;
  double beta_now = 0;
  
  if(v_now<v_out){
    beta_now = beta_out;
  }else{
    beta_now = v_now/300.;
  }

  // std::cout << "recoil a " << a << " t " << (v_out-v_in)/a << "ns " << v_now << " " << beta_now << std::endl;
  return beta_now;
}

//
double angle_info::ComputeCosThetaParticleGamma(double tg, double pg, double tp, double pp)
{
  tg = tg/180.*PI;
  pg = pg/180.*PI;
  tp = tp/180.*PI;
  pp = pp/180.*PI;

  double tmp1 = std::sin(tp)*std::sin(tg)*std::cos(pp-pg);
  double tmp2 = std::cos(tp)*std::cos(tg);

  return tmp1+tmp2;
}

namespace sm150_doppler {
namespace {
constexpr double rad = 3.14159265358979323846 / 180.0;
// Preserve the historical mass convention for a controlled baseline comparison.
constexpr double mass_p = 16*938.3 + 16*931.5 - 32*8.493;
constexpr double mass_r = 62*938.3 + 88*931.5 - 150*8.262;
void check(bool condition,const std::string &message) {
  if(!condition) throw std::runtime_error(message);
}
double cosine(double tg,double pg,double tr,double pr) {
  return std::sin(tg*rad)*std::sin(tr*rad)*std::cos((pg-pr)*rad)+std::cos(tg*rad)*std::cos(tr*rad);
}
}
Model::Model(const std::string &project,const nlohmann::json &parameters):config(parameters) {
  const sm150::ParameterPaths runtime(std::filesystem::path(project)/"pars");
  const int schema=config.at("schema_version");
  check(schema>=1&&schema<=3&&config.at("schema_version")==schema,"Unsupported Doppler parameter schema");
  tabulated_beta=schema==3;
  check(schema==2||!config.contains("energy_scale_by_period"),"Period energy scales require Doppler schema2");
  check(tabulated_beta||!config.contains("beta_by_period"),"Direct beta tables require Doppler schema3");
  for(const auto &family:{"Si","S3"}) {
    auto &x=config.at(family);
    const double theta=x.at("theta_offset_deg"),phi=x.at("phi_offset_deg"),blend=x.at("particle_energy_fraction");
    if(tabulated_beta) {
      check(!x.contains("energy_scale"),"Direct beta must not also apply an energy scale");
      check(theta==0&&blend==0,"Ring beta tables require nominal particle theta and no event-energy blend");
    } else {
      const double scale=x.at("energy_scale");
      check(std::isfinite(scale)&&scale>0&&scale<=3,"Invalid recoil energy scale");
    }
    check(std::isfinite(theta)&&std::abs(theta)<=5&&std::isfinite(phi)&&std::abs(phi)<=30,"Geometry scan exceeds allowed range");
    check(std::isfinite(blend)&&blend>=0&&blend<=1,"Invalid particle-energy fraction");
  }
  // Schema1 keeps one scale for all periods. Schema2 can override only the
  // requested family and period; old executables reject it instead of silently
  // treating a period-dependent configuration as a global baseline.
  if(!tabulated_beta)for(auto &scales:energy_scales)scales={config.at("Si").at("energy_scale"),config.at("S3").at("energy_scale")};
  if(config.contains("energy_scale_by_period")) {
    const auto &overrides=config.at("energy_scale_by_period");
    check(overrides.is_object(),"Period energy scales must be an object");
    const std::array<std::string,3> periods={"high1","low","high2"};
    for(auto it=overrides.begin();it!=overrides.end();++it) {
      const auto where=std::find(periods.begin(),periods.end(),it.key());
      check(where!=periods.end()&&it.value().is_object(),"Invalid Doppler period override");
      const auto group=std::size_t(where-periods.begin());
      for(auto family=it.value().begin();family!=it.value().end();++family) {
        check(family.key()=="Si"||family.key()=="S3","Unknown particle family in period energy scales");
        const double scale=family.value();
        check(std::isfinite(scale)&&scale>0&&scale<=3,"Invalid period recoil energy scale");
        energy_scales[group][family.key()=="S3" ? 1:0]=scale;
      }
    }
  }
  auto open=[&](const std::string &name) {
    std::ifstream f(runtime.legacy(name));
    check(bool(f),"Cannot read "+name); return f;
  };
  for(int which=0;which<2;++which) {
    auto f=open(which ? "angle_info/spider_angle.txt":"angle_info/ge_angle.txt");
    std::string header;std::getline(f,header);int r,s;double t,p;
    auto &out=which ? si_angles:ge_angles;
    while(f>>r>>s>>t>>p) check(out.emplace(100*r+s,std::make_pair(t,p)).second,"Duplicate angle entry");
    check(f.eof()&&!out.empty(),"Malformed detector angle table");
  }
  for(int which=0;which<2;++which) {
    auto f=open(which ? "angle_info/s3_sector_angle.txt":"angle_info/s3_ring_angle.txt");
    std::string header;std::getline(f,header);int id;double angle;
    auto &out=which ? s3_phi:s3_theta;
    while(f>>id>>angle) check(out.emplace(id,angle).second,"Duplicate S3 angle entry");
    check(f.eof()&&!out.empty(),"Malformed S3 angle table");
  }
  if(tabulated_beta) {
    check(config.at("beta_units")=="v/c","Unexpected beta units");
    const auto &hashes=config.at("beta_input_sha256");
    check(hashes.is_object()&&hashes.size()==4,"Incomplete beta input provenance");
    for(const auto *name:{"angle_info/spider_angle.txt","angle_info/s3_ring_angle.txt", "beta_info/theta_ein_eout_p.txt","beta_info/theta_ein_eout_r.txt"})
      check(hashes.at(name)==sm150::ParameterPaths::sha256(runtime.legacy(name)),"Beta input SHA256 mismatch: "+std::string(name));
    std::array<std::map<int,double>,2> polar_angles;
    polar_angles[1]=s3_theta;
    for(const auto &item:si_angles) {
      auto inserted=polar_angles[0].emplace(item.first/100,item.second.first);
      check(inserted.second||inserted.first->second==item.second.first,"Si ring has sector-dependent polar angles");
    }
    const auto &tables=config.at("beta_by_period");
    check(tables.is_object()&&tables.size()==3,"Require all three beta periods");
    const std::array<std::string,3> periods={"high1","low","high2"};
    for(int group=0;group<3;++group) {
      const auto &period=tables.at(periods[group]);
      check(period.is_object()&&period.size()==2,"Require both beta families");
      for(int ref=0;ref<2;++ref) {
        const auto &rings=period.at(ref ? "S3":"Si");
        check(rings.is_object()&&rings.size()==polar_angles[ref].size(),"Incomplete ring beta table");
        for(const auto &ring:polar_angles[ref]) {
          const auto &entry=rings.at(std::to_string(ring.first));
          RingKinematics k{entry.at("beta"),entry.at("theta_p_deg"),entry.at("theta_r_deg"),
            entry.at("expected_particle_MeV"),entry.at("kinematic_extrapolated")};
          check(std::isfinite(k.beta)&&k.beta>=0&&k.beta<1,"Invalid tabulated beta");
          check(k.theta_p==ring.second,"Beta table and particle geometry disagree");
          check(std::isfinite(k.theta_r)&&k.theta_r>=0&&k.theta_r<=180&&
            std::isfinite(k.expected_particle)&&k.expected_particle>0,"Invalid ring kinematics");
          ring_kinematics[group][ref].emplace(ring.first,k);
        }
      }
    }
    return;
  }
  std::map<int,std::pair<double,double>> projectile;
  auto fp=open("beta_info/theta_ein_eout_p.txt");double p,r,ei,eo;
  while(fp>>p>>ei>>eo) check(projectile.emplace(int(std::lround(p*100)),std::make_pair(ei,eo)).second,"Duplicate projectile angle");
  check(fp.eof(),"Malformed projectile table");
  auto fr=open("beta_info/theta_ein_eout_r.txt");
  while(fr>>p>>r>>ei>>eo) {
    auto e=projectile.at(int(std::lround(p*100)));
    knots[p<90 ? 1:0].push_back({std::round(p*100)/100,std::round(r*100)/100,e.first,e.second,ei,eo});
  }
  check(fr.eof(),"Malformed recoil table");
  for(auto &table:knots) {
    std::sort(table.begin(),table.end(),[](const Knot &a,const Knot &b){return a.p<b.p;});
    check(table.size()>=2,"Missing kinematic table");
  }
}
double Model::energy_scale(int group,int reference) const {
  check(!tabulated_beta,"Direct beta configuration has no runtime energy scale");
  check(group>=0&&group<3,"Invalid Doppler data period");
  check(reference==0||reference==1,"Invalid particle family");
  return energy_scales.at(group).at(reference);
}
Model::Prepared Model::prepare(const Pair &pair) {
  check(pair.reference==0||pair.reference==1,"Invalid particle family");
  check(pair.group>=0&&pair.group<3,"Invalid Doppler data period");
  const auto &par=config.at(pair.reference ? "S3":"Si");
  auto gamma=ge_angles.at(100*pair.ge_ring+pair.ge_sector);
  auto particle=pair.reference ? std::make_pair(s3_theta.at(pair.ring),s3_phi.at(pair.sector)) : si_angles.at(100*pair.ring+pair.sector);
  const double nominal=std::round(particle.first*100)/100;
  particle.first+=par.at("theta_offset_deg").get<double>();
  particle.second+=par.at("phi_offset_deg").get<double>();
  gamma.first+=config.value("ge_theta_offset_deg",0.0);
  gamma.second+=config.value("ge_phi_offset_deg",0.0);
  if(config.contains("per_ge")&&config["per_ge"].contains(std::to_string(pair.ge))) {
    const auto &g=config["per_ge"][std::to_string(pair.ge)];
    gamma.first+=g.value("theta_offset_deg",0.0);gamma.second+=g.value("phi_offset_deg",0.0);
  }
  Prepared z{};
  z.c.theta_p=particle.first;z.c.phi_p=particle.second;
  z.c.phi_r=particle.second+180;z.c.theta_g=gamma.first;z.c.phi_g=gamma.second;
  if(tabulated_beta) {
    const auto &k=ring_kinematics.at(pair.group).at(pair.reference).at(pair.ring);
    z.c.beta=k.beta;z.c.theta_r=k.theta_r;z.c.expected_particle_MeV=k.expected_particle;
    z.c.extrapolated=k.extrapolated;
    z.c.cosine=cosine(gamma.first,gamma.second,k.theta_r,z.c.phi_r);
    return z;
  }
  const auto &table=knots.at(pair.reference);
  const double x=nominal+par.at("theta_offset_deg").get<double>();
  auto upper=std::upper_bound(table.begin(),table.end(),x,[](double t,const Knot &k){return t<k.p;});
  size_t j=upper==table.begin() ? 1 : std::min(size_t(upper-table.begin()),table.size()-1);
  const auto &a=table[j-1],&b=table[j];const double f=(x-a.p)/(b.p-a.p);
  auto mix=[&](double u,double v){return u+f*(v-u);};
  const double tr=mix(a.r,b.r),pi=mix(a.pi,b.pi),po=mix(a.po,b.po),ri=mix(a.ri,b.ri),ro=mix(a.ro,b.ro);
  check(pi>0&&po>0&&ri>0&&ro>0,"Invalid interpolated kinematics");
  z.c.theta_r=tr;z.c.expected_particle_MeV=po;
  z.c.extrapolated=x<table.front().p||x>table.back().p;
  z.scale=energy_scale(pair.group,pair.reference);z.energy_out=ro;
  z.blend=par.at("particle_energy_fraction");
  z.pp_base=std::sqrt(2*mass_p*pi);z.beam_p=std::sqrt(2*mass_p*(pi+ri));
  const double px=z.pp_base*std::sin(particle.first*rad),pz=z.beam_p-z.pp_base*std::cos(particle.first*rad);
  z.pr_base=(px*px+pz*pz)/(2*mass_r);z.pt_base=std::atan2(px,pz)/rad;
  z.c.beta=std::sqrt(2*ro*z.scale/mass_r);
  z.c.cosine=cosine(gamma.first,gamma.second,tr,z.c.phi_r);
  return z;
}
Correction Model::evaluate(const Pair &pair) {
  const auto key=std::make_tuple(pair.group,pair.ge,pair.ge_ring,pair.ge_sector,pair.reference,pair.ring,pair.sector);
  auto it=cache.find(key);if(it==cache.end()) it=cache.emplace(key,prepare(pair)).first;
  const auto &z=it->second;Correction c=z.c;
  if(z.blend>0) {
    // A table-anchored momentum perturbation, not an independently calibrated velocity.
    const double ep=z.pp_base*z.pp_base/(2*mass_p)+z.blend*(pair.particle_MeV-c.expected_particle_MeV);
    check(ep>0,"Nonpositive reconstructed projectile energy");
    const double pp=std::sqrt(2*mass_p*ep),px=pp*std::sin(c.theta_p*rad),pz=z.beam_p-pp*std::cos(c.theta_p*rad);
    const double exit=z.energy_out+(px*px+pz*pz)/(2*mass_r)-z.pr_base;
    c.stopped=exit<0;c.beta=std::sqrt(2*std::max(0.0,exit)*z.scale/mass_r);
    c.theta_r+=std::atan2(px,pz)/rad-z.pt_base;
    c.cosine=cosine(c.theta_g,c.phi_g,c.theta_r,c.phi_r);
  }
  check(std::isfinite(c.beta)&&c.beta>=0&&c.beta<1&&std::abs(c.cosine)<=1.0000001,"Invalid Doppler factor");
  c.energy=pair.raw*(1-c.beta*c.cosine)/std::sqrt(1-c.beta*c.beta);
  return c;
}
}

