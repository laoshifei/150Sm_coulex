#ifndef _doppler_h_
#define _doppler_h_

#include "set.h"
#include "event.h"
#include "angle_info.h"

#include <vector>
#include <map>

#include "TH1D.h"
#include "TBenchmark.h"

//
class doppler
{
public:
  doppler(const std::string &filename_in, const std::string &filename_out, std::string config);
  virtual ~doppler();

  void Process();
  void ProcessDoppler(event *ed, TTree *tr);
  void SaveFile();

public:
  double ComputeEnergyWithDoppler(double adc_ge, double beta, double tg, double pg, double tp, double pp);
  double ComputeEnergyWithDoppler(double adc_ge, double beta, double cospg);

private:
  bool ReadCutInfo();
  void PrintCutInfo();

  bool ReadGeConfigInfo();
  void PrintGeConfigInfo();

private:
  void Clear();

private:
  angle_info *angle;

private:
  std::map<int, bool> map_ge_flag;

  std::map<int, std::pair<double, double>> map_spider_cut_p;
  std::map<int, std::pair<double, double>> map_spider_cut_r;

  std::map<int, std::pair<double, double>> map_s3_sector_cut_p;
  std::map<int, std::pair<double, double>> map_s3_sector_cut_r;
  std::map<int, std::pair<double, double>> map_s3_ring_cut_p;
  std::map<int, std::pair<double, double>> map_s3_ring_cut_r;

private:
  std::string configure;

private:
  TBenchmark *benchmark;
  TFile *file_in;
  TTree *tr_event;
  TTree *tr_bg;

  event *event_coin;
  event *event_bg;

  TFile *file_out;
  TTree *tr_event_out;
  TTree *tr_bg_out;
  
  Int_t hits;
  Bool_t flag_spider_p;
  Bool_t flag_spider_r;
  Bool_t flag_s3_p;
  Bool_t flag_s3_r;
  Int_t rp; // ring for projectile
  Int_t sp; // sector for projectile
  Int_t rr; // ring for recoil
  Int_t sr; // sector for recoil
  Int_t rg[GENUM]; // ring for gamma
  Int_t sg[GENUM]; // sector for gamma

  Double_t e_raw[GENUM];
  Double_t e_dc_p[GENUM];
  Double_t e_dc_r[GENUM];
};

namespace sm150_doppler {
struct Options {
  std::string project,parameters,input,output,scan_list,scan_parameters;
};
int process_build(const Options &options);
int scan_pairs(const Options &options);
}
#endif
