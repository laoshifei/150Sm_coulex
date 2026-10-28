#ifndef _build_h_
#define _build_h_

#include "set.h"
#include "../include/sm150_time_alignment.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TBenchmark.h"
#include "TSystem.h"

#include <iostream>

class build
{
public:
  build(const std::string &filename_in, const std::string &filename_out, int r);
  virtual ~build();

  void Process();

  void GetGeSpiderS3EventPrompt();
  void GetGeSpiderS3EventRandom();
  void SaveFile();

private:
  std::map<int, std::vector<double>> map_cali_data;
  std::map<int, double> map_ts_offset;
  std::map<int, std::vector<double>> map_ge_tw_data;
  sm150::TimeAlignment time_alignment;
  double timing_ge_energy_keV = 0.;

  std::map<int, double> map_k_data;

  std::map<int, int> map_s3_sector_id, map_s3_ring_id;
  std::map<int, std::vector<double>> map_s3_sector_cor_data, map_s3_ring_cor_data;

private:
  bool ReadCaliData();
  void PrintCaliData();

  void ReadS3CorData();
  void PrintS3CorData();

  bool ReadTsOffset();
  void PrintTsOffset();

  bool ReadGeTwData();
  void PrintGeTwData();

  bool InitMapSectorRingID();
  void PrintMapSectorRingID();

private:
  double GetEnergy();
  Long64_t GetTsns();

private:
  TBenchmark *benchmark;  
  TRandom3 *rndm;  

private:
  Int_t run;
  std::vector<int> v_s3_sid;

private:
  TFile *file_in;
  TTree *tr_in;
  Short_t sr;
  Short_t cid;
  Short_t sid;
  Short_t ch;
  UShort_t evte;
  Long64_t ts;
  Short_t cfd;
  Bool_t cfdft;
  Short_t cfds;

  Long64_t ts_ns;
  Double_t energy;

  TF1 *tf;

  TFile *file_out;
};

#endif
