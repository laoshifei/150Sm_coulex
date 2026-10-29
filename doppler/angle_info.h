#ifndef angle_info_h
#define angle_info_h

#define PI 3.14159265

#include <map>
#include <vector>
#include <string>
#include <array>
#include <tuple>
#include <nlohmann/json.hpp>

#include "TString.h"

//
class angle_info{
public:
  angle_info();
  virtual ~angle_info();

private:
  std::map<int, std::pair<double, double>> map_ge_angle; // theta,phi
  std::map<int, std::pair<double, double>> map_spider_angle; // theta,phi
  std::map<int, double> map_s3_sector_phi;
  std::map<int, double> map_s3_ring_theta;

  std::map<std::string, std::pair<double, double>> map_theta_ein_eout_p;
  std::map<std::string, std::pair<double, double>> map_theta_ein_eout_r;
  std::map<std::string, double> map_theta_p2r;

public:
  bool InitAngleInfo();
  void PrintAngleInfo();

  bool InitEinEoutInfo();
  void PrintEinEoutInfo();
  void PrintThetaP2RInfo();

  void PrintBetaInfo();

public:
  double GetGammaTheta(int ring, int sector) { return map_ge_angle[100*ring+sector].first; }
  double GetGammaPhi(int ring, int sector) { return map_ge_angle[100*ring+sector].second; }
  
  double GetProjectileSpiderTheta(int ring, int sector) { return map_spider_angle[100*ring+sector].first; } 
  double GetProjectileSpiderPhi(int ring, int sector) { return map_spider_angle[100*ring+sector].second; }
  double GetProjectileS3Theta(int ring) { return map_s3_ring_theta[ring]; } 
  double GetProjectileS3Phi(int sector) { return map_s3_sector_phi[sector]; }

  double GetRecoilTheta(double theta_p);
  double GetRecoilPhi(double phi_p);

  double ComputeProjectileBeta(double theta, double dt);
  double ComputeRecoilBeta(double theta, double dt); // theta

public:
  // tg: theta gamma
  // pg: phi gamma
  // tp: theta particle
  // pp: phi particle
  double ComputeCosThetaParticleGamma(double tg, double pg, double tp, double pp);
};

namespace sm150_doppler {
// One accepted Ge-particle pair; detector energies are already calibrated.
struct Pair {
  int run=0, group=0, ge=0, ge_ring=0, ge_sector=0, reference=0;
  int particle=0, ring=0, sector=0, side=0, ge_flags=0;
  long long ge_source=0, particle_source=0, sector_source=-1;
  double raw=0, particle_MeV=0, dt_ns=0, alpha=0, extra_factor=1;
  bool legacy_selection=false;
};
struct Correction {
  double energy=0, beta=0, cosine=0, theta_p=0, phi_p=0, theta_r=0, phi_r=0;
  double theta_g=0, phi_g=0, expected_particle_MeV=0;
  bool extrapolated=false, stopped=false;
};
class Model {
public:
  Model(const std::string &project, const nlohmann::json &parameters);
  Correction evaluate(const Pair &pair);
  const nlohmann::json &parameters() const { return config; }
private:
  struct Knot { double p,r,pi,po,ri,ro; };
  struct Prepared { Correction c; double energy_out,pr_base,pt_base,beam_p,pp_base,blend,scale; };
  nlohmann::json config;
  std::map<int,std::pair<double,double>> ge_angles,si_angles;
  std::map<int,double> s3_theta,s3_phi;
  std::array<std::vector<Knot>,2> knots;
  std::map<std::tuple<int,int,int,int,int,int>,Prepared> cache;
  Prepared prepare(const Pair &pair);
};
}

#endif
