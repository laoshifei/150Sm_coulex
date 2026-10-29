#ifndef SORT_H_
#define SORT_H_

#include "set.h"

#include "TBenchmark.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TRandom3.h"

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <ctime>

#include "raw.h"
#include "../include/sm150_time_alignment.h"

//
struct sort_data {
  Short_t cid;
  Short_t sid;
  Short_t ch;
  Double_t evte;
  Long64_t ts;

  sort_data() : cid(0), sid(0), ch(0), evte(0.), ts(0) {}

  void operator()(Short_t cid_, Short_t sid_, Short_t ch_, Double_t evte_, Long64_t ts_)
  {
    cid = cid_;
    sid = sid_;
    ch = ch_;
    evte = evte_;
    ts = ts_;
  }
};

//
class sort
{
public:
  sort() = default;
  sort(const std::string &filename_in, const std::string &filename_out, const int r);
  ~sort();

public:
  bool ReadCaliData();
  void PrintCaliData();

  bool ReadS3CorData();
  void PrintS3CorData();

  bool ReadS3AttData();
  void PrintS3AttData();

  bool ReadTSOffset();
  void PrintTSOffset();

  bool ReadGeTwData();
  void PrintGeTwData();

  bool InitMapSectorRingID();
  void PrintMapSectorRingID();

public:
  void Process();

private:
  void GetEnergy();
  void GetTsns();

private:
  std::map<int, std::vector<double>> map_cali;
  std::map<int, double> map_k;
  std::map<int, std::vector<double>> map_s3_sector_cor, map_s3_ring_cor; // key is id
  std::map<int, double> map_s3_sector_att, map_s3_ring_att; // key is id
  std::map<int, double> map_ts_offset;
  std::map<int, std::vector<double>> map_ge_tw;
  sm150::TimeAlignment time_alignment;
  double timing_ge_energy_keV = 0.;

  std::map<int, int> map_s3_sector_id, map_s3_ring_id;

  std::map<Long64_t, sort_data> map_sort_data;
  std::map<Long64_t, sort_data>::iterator it_sd;

private:
  TBenchmark *benchmark;

  int run;
  TFile *file_in;
  TTree *tr_in;

  raw *rd;

  TF1 *tf;
  TRandom3 *rndm;

private:
  TFile *file_out;

  TTree *tr_out;
  sort_data sd;
};

#endif
