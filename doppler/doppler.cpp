#include "doppler.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include "TH1D.h"

#include "set.h"
#include <TObjString.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <openssl/evp.h>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <array>
#include <cmath>

//
doppler::doppler(const std::string &filename_in, const std::string &filename_out, std::string config)
{
  configure = config;

  benchmark = new TBenchmark;

  file_in = TFile::Open(filename_in.c_str());
  if(file_in->IsZombie()){
    std::cout << "open file " << filename_in.c_str() << " error!" << std::endl;
    delete file_in;

    throw std::invalid_argument("can not open rootfile to read.");
  }

  file_out = new TFile(filename_out.c_str(), "recreate");
  if(!file_out){
    std::cout << "can not create " << filename_out << std::endl;
    delete file_out;
    
    throw std::invalid_argument("can not open rootfile to write.");
  }
  
  tr_event_out = new TTree("tr_event", "event doppler data");
  tr_event_out->Branch("hits", &hits, "hits/I");
  tr_event_out->Branch("flag_spider_p", &flag_spider_p, "flag_spider_p/O");
  tr_event_out->Branch("flag_spider_r", &flag_spider_r, "flag_spider_r/O");
  tr_event_out->Branch("flag_s3_p", &flag_s3_p, "flag_s3_p/O");
  tr_event_out->Branch("flag_s3_r", &flag_s3_r, "flag_s3_r/O");
  tr_event_out->Branch("rp", &rp, "rp/I");
  tr_event_out->Branch("sp", &sp, "sp/I");
  tr_event_out->Branch("rr", &rr, "rr/I");
  tr_event_out->Branch("sr", &sr, "sr/I");
  tr_event_out->Branch("rg", rg, "rg[hits]/I");
  tr_event_out->Branch("sg", sg, "sg[hits]/I");

  tr_event_out->Branch("e_raw", e_raw, "e_raw[hits]/D");
  tr_event_out->Branch("e_dc_p", e_dc_p, "e_dc_p[hits]/D");
  tr_event_out->Branch("e_dc_r", e_dc_r, "e_dc_r[hits]/D");

  tr_bg_out = new TTree("tr_bg", "bg doppler data");
  tr_bg_out->Branch("hits", &hits, "hits/I");
  tr_bg_out->Branch("flag_spider_p", &flag_spider_p, "flag_spider_p/O");
  tr_bg_out->Branch("flag_spider_r", &flag_spider_r, "flag_spider_r/O");
  tr_bg_out->Branch("flag_s3_p", &flag_s3_p, "flag_s3_p/O");
  tr_bg_out->Branch("flag_s3_r", &flag_s3_r, "flag_s3_r/O");
  tr_bg_out->Branch("rp", &rp, "rp/I");
  tr_bg_out->Branch("sp", &sp, "sp/I");
  tr_bg_out->Branch("rr", &rr, "rr/I");
  tr_bg_out->Branch("sr", &sr, "sr/I");
  tr_bg_out->Branch("rg", rg, "rg[hits]/I");
  tr_bg_out->Branch("sg", sg, "sg[hits]/I");

  tr_bg_out->Branch("e_raw", e_raw, "e_raw[hits]/D");
  tr_bg_out->Branch("e_dc_p", e_dc_p, "e_dc_p[hits]/D");
  tr_bg_out->Branch("e_dc_r", e_dc_r, "e_dc_r[hits]/D");

  //
  tr_event = (TTree*)file_in->Get("tr_event");
  tr_bg = (TTree*)file_in->Get("tr_bg");
  event_coin = new event(tr_event); 
  event_bg = new event(tr_bg); 
  
  //angle info
  angle = new angle_info();

  //
  if(!ReadGeConfigInfo()){
    throw std::invalid_argument("can not read ge configure info data.");
  }
  PrintGeConfigInfo();

  // cut info
  if(!ReadCutInfo()){
    throw std::invalid_argument("can not read cut info data.");
  }
  PrintCutInfo();
}

//
doppler::~doppler()
{
  delete file_in;
  delete event_coin;
  delete event_bg;
  delete file_out;
}

//
void doppler::Process()
{
  benchmark->Start("build");

  ProcessDoppler(event_coin, tr_event_out);
  ProcessDoppler(event_bg, tr_bg_out);

  benchmark->Show("doppler");
}

//
void doppler::ProcessDoppler(event *ed, TTree *tr)
{
  file_out->cd();
  Clear();

  std::vector<int> v_ge_index;
  std::vector<int> v_spider_p_index;
  std::vector<int> v_spider_r_index;
  std::vector<int> v_s3_sector_p_index;
  std::vector<int> v_s3_sector_r_index;
  std::vector<int> v_s3_ring_p_index;
  std::vector<int> v_s3_ring_r_index;

  //
  int num_spider_p = 0, num_spider_r = 0;
  int num_s3_sector_p = 0, num_s3_sector_r = 0;
  int num_s3_ring_p = 0, num_s3_ring_r = 0;

  int key = 0;

  double tg, pg, tp_p, pp_p, tp_r, pp_r;
  double cospg_p, cospg_r;
  double beta_p, beta_r;

  std::cout << "\ntotal entries " << ed->GetEntries() << std::endl;
  for(Long64_t i=0;i<ed->GetEntries();++i){
    if(i%10000==0){ 
      std::cout << "\r" << i << "/" << ed->GetEntries();
      std::cout << std::flush;
    }

    Clear();

    v_ge_index.clear();
    v_spider_p_index.clear();
    v_spider_r_index.clear();
    v_s3_sector_p_index.clear();
    v_s3_sector_r_index.clear();
    v_s3_ring_p_index.clear();
    v_s3_ring_r_index.clear();

    num_spider_p = 0; num_spider_r = 0; 
    num_s3_sector_p = 0; num_s3_sector_r = 0;
    num_s3_ring_p = 0; num_s3_ring_r = 0;
    flag_spider_p = 0; flag_spider_r = 0;
    flag_s3_p = 0; flag_s3_r = 0;

    key = 0;
    ed->GetEntry(i);

    // ge
    for(int j=0;j<ed->n_ge;j++){
      if(map_ge_flag[100*ed->ge_sid[j]+ed->ge_ch[j]]){
        v_ge_index.push_back(j);
      }
    }

    // spider
    for(int j=0;j<ed->n_spider;j++){
      key = 100*ed->spider_ring_id[j]+ed->spider_sector_id[j];
      if(ed->spider_energy[j]>map_spider_cut_p[key].first && ed->spider_energy[j]<map_spider_cut_p[key].second){
        num_spider_p++;
        v_spider_p_index.push_back(j);
      }
      if(ed->spider_energy[j]>map_spider_cut_r[key].first && ed->spider_energy[j]<map_spider_cut_r[key].second){
        num_spider_r++;
        v_spider_r_index.push_back(j);
      }
    }

    // s3 sector
    for(int j=0;j<ed->n_s3_sector;j++){
      key = ed->s3_sector_id[j];
      if(ed->s3_sector_energy[j]>map_s3_sector_cut_p[key].first && ed->s3_sector_energy[j]<map_s3_sector_cut_p[key].second){
        num_s3_sector_p++;
        v_s3_sector_p_index.push_back(j);
      }
      if(ed->s3_sector_energy[j]>map_s3_sector_cut_r[key].first && ed->s3_sector_energy[j]<map_s3_sector_cut_r[key].second){
        num_s3_sector_r++;
        v_s3_sector_r_index.push_back(j);
      }
    }

    // s3 ring
    for(int j=0;j<ed->n_s3_ring;j++){
      key = ed->s3_ring_id[j];
      if(ed->s3_ring_energy[j]>map_s3_ring_cut_p[key].first && ed->s3_ring_energy[j]<map_s3_ring_cut_p[key].second){
        num_s3_ring_p++;
        v_s3_ring_p_index.push_back(j);
      }
      if(ed->s3_ring_energy[j]>map_s3_ring_cut_r[key].first && ed->s3_ring_energy[j]<map_s3_ring_cut_r[key].second){
        num_s3_ring_r++;
        v_s3_ring_r_index.push_back(j);
      }
    }

    // cout
    // std::cout << "\ni " << i 
    //           << " \nge enrgy             " << ed->ge_energy[0]
    //           << " \nnumber_spider_p      " << num_spider_p
    //           << " \nnumber_spider_r      " << num_spider_r
    //           << " \nnumber_s3_sector_p   " << num_s3_sector_p
    //           << " \nnumber_s3_sector_r   " << num_s3_sector_r
    //           << " \nnumber_s3_ring_p     " << num_s3_ring_p
    //           << " \nnumber_s3_ring_r     " << num_s3_ring_r
    //           << " \nflag_ge size          " << v_ge_index.size()
    //           << " \nflag_spider_p size    " << v_spider_p_index.size()
    //           << " flag_spider_r size      " << v_spider_r_index.size()
    //           << " \nflag_s3_sector_p size " << v_s3_sector_p_index.size()
    //           << " flag_s3_sector_r size   " << v_s3_sector_r_index.size()
    //           << " \nflag_s3_ring_p size   " << v_s3_ring_p_index.size()
    //           << " flag_s3_ring_r size     " << v_s3_ring_r_index.size()
    //           << "\n" << std::endl;
    // for(int j=0;j<ed->n_ge;j++){
    //   std::cout << "ge energy " << ed->ge_energy[j] << std::endl;
    //   std::cout << "ge ring id " << ed->ge_ring_id[j] << std::endl;
    //   std::cout << "ge sector id " << ed->ge_sector_id[j] << std::endl;
    // }
    // for(int j=0;j<ed->n_spider;j++){
    //   std::cout << "spider energy " << ed->spider_energy[j] << std::endl;
    //   std::cout << "spider ring id " << ed->spider_ring_id[j] << std::endl;
    //   std::cout << "spider sector id " << ed->spider_sector_id[j] << std::endl;
    // }
    // for(int j=0;j<ed->n_s3_ring;j++){
    //   std::cout << "s3 ring energy " << ed->s3_ring_energy[j] << std::endl;
    //   std::cout << "s3 ring id " << ed->s3_ring_id[j] << std::endl;
    // }
    // for(int j=0;j<ed->n_s3_sector;j++){
    //   std::cout << "s3 sector energy " << ed->s3_sector_energy[j] << std::endl;
    //   std::cout << "s3 sector id " << ed->s3_sector_id[j] << std::endl;
    // }
    // std::cout << "\n\n";

    // spider p
    if(num_spider_p>1){
      int max_index = v_spider_p_index[0];
      double max_e = ed->spider_energy[max_index];

      for(int j=1;j<num_spider_p;j++){
        int idx = v_spider_p_index[j];
        if(ed->spider_energy[idx]>max_e){
          max_index = idx;
          max_e = ed->spider_energy[idx];
        }
      }

      v_spider_p_index.clear();
      v_spider_p_index.push_back(max_index);
      num_spider_p = 1;
    }
    if(num_spider_p==1) flag_spider_p = 1;

    // spider r
    if(num_spider_r>1){
      int max_index = v_spider_r_index[0];
      double max_e = ed->spider_energy[max_index];

      for(int j=1;j<num_spider_r;j++){
        int idx = v_spider_r_index[j];
        if(ed->spider_energy[idx]>max_e){
          max_index = idx;
          max_e = ed->spider_energy[idx];
        }
      }

      v_spider_r_index.clear();
      v_spider_r_index.push_back(max_index);
      num_spider_r = 1;
    }
    if(num_spider_r==1) flag_spider_r = 1;

    // s3 p
    if(num_s3_sector_p==0 || num_s3_ring_p==0){
      flag_s3_p = 0;
      num_s3_sector_p = 0;
      num_s3_ring_p = 0;
    }else if (num_s3_sector_p==1 && num_s3_ring_p==1){
      flag_s3_p = 1;
    }else if (num_s3_sector_p==1 && num_s3_ring_p==2){
      flag_s3_p = 1;
      num_s3_sector_p = 1;
      num_s3_ring_p = 1;

      double e1 = ed->s3_sector_energy[v_s3_sector_p_index[0]];
      double e2 = ed->s3_ring_energy[v_s3_ring_p_index[0]];
      double e3 = ed->s3_ring_energy[v_s3_ring_p_index[1]];

      if(std::abs(e1-e2) < (std::abs(e1-e3))){
        v_s3_ring_p_index.erase(v_s3_ring_p_index.begin()+1);
      }else{
        v_s3_ring_p_index.erase(v_s3_ring_p_index.begin());
      }
    }else if(num_s3_sector_p==2 && num_s3_ring_p==1){
      flag_s3_p = 1;
      num_s3_sector_p = 1;
      num_s3_ring_p = 1;

      double e1 = ed->s3_ring_energy[v_s3_ring_p_index[0]];
      double e2 = ed->s3_sector_energy[v_s3_sector_p_index[0]];
      double e3 = ed->s3_sector_energy[v_s3_sector_p_index[1]];

      if(std::abs(e1-e2) < (std::abs(e1-e3))){
        v_s3_sector_p_index.erase(v_s3_sector_p_index.begin()+1);
      }else{
        v_s3_sector_p_index.erase(v_s3_sector_p_index.begin());
      }
    }else{
      flag_s3_p = 0;
      num_s3_sector_p = 0;
      num_s3_ring_p = 0;

      v_s3_sector_p_index.clear();
      v_s3_ring_p_index.clear();
    }

    // s3 r
    if(num_s3_sector_r==0 || num_s3_ring_r==0){
      flag_s3_r = 0;
      num_s3_sector_r = 0;
      num_s3_ring_r = 0;
    }else if (num_s3_sector_r==1 && num_s3_ring_r==1){
      flag_s3_r = 1;
    }else if (num_s3_sector_r==1 && num_s3_ring_r==2){
      flag_s3_r = 1;
      num_s3_sector_r = 1;
      num_s3_ring_r = 1;

      double e1 = ed->s3_sector_energy[v_s3_sector_r_index[0]];
      double e2 = ed->s3_ring_energy[v_s3_ring_r_index[0]];
      double e3 = ed->s3_ring_energy[v_s3_ring_r_index[1]];

      if(std::abs(e1-e2) < (std::abs(e1-e3))){
        v_s3_ring_r_index.erase(v_s3_ring_r_index.begin()+1);
      }else{
        v_s3_ring_r_index.erase(v_s3_ring_r_index.begin());
      }
    }else if(num_s3_sector_r==2 && num_s3_ring_r==1){
      flag_s3_r = 1;
      num_s3_sector_r = 1;
      num_s3_ring_r = 1;

      double e1 = ed->s3_ring_energy[v_s3_ring_r_index[0]];
      double e2 = ed->s3_sector_energy[v_s3_sector_r_index[0]];
      double e3 = ed->s3_sector_energy[v_s3_sector_r_index[1]];

      if(std::abs(e1-e2) < (std::abs(e1-e3))){
        v_s3_sector_r_index.erase(v_s3_sector_r_index.begin()+1);
      }else{
        v_s3_sector_r_index.erase(v_s3_sector_r_index.begin());
      }
    }else{
      flag_s3_r = 0;
      num_s3_sector_r = 0;
      num_s3_ring_r = 0;

      v_s3_sector_r_index.clear();
      v_s3_ring_r_index.clear();
    }

    // cout
    // for(int j=0;j<ed->n_ge;j++){
    //   std::cout << "ge energy " << ed->ge_energy[j] << std::endl;
    //   std::cout << "ge theta " << angle->GetGammaTheta(ed->ge_ring_id[j], ed->ge_sector_id[j]) << std::endl;
    //   std::cout << "ge phi " << angle->GetGammaPhi(ed->ge_ring_id[j], ed->ge_sector_id[j]) << std::endl;
    // }
    // if(flag_spider_p){
    //   for(int j=0;j<num_spider_p;j++){
    //     std::cout << "spider p energy " << ed->spider_energy[v_spider_p_index[j]] << std::endl;
    //     std::cout << "spider p theta " << angle->GetProjectileSpiderTheta(ed->spider_ring_id[v_spider_p_index[j]], ed->spider_sector_id[v_spider_p_index[j]]) << std::endl;
    //     std::cout << "spider p phi " << angle->GetProjectileSpiderPhi(ed->spider_ring_id[v_spider_p_index[j]], ed->spider_sector_id[v_spider_p_index[j]]) << std::endl;
    //   }
    // }
    // if(flag_spider_r){
    //   for(int j=0;j<num_spider_r;j++){
    //     std::cout << "spider r energy " << ed->spider_energy[v_spider_r_index[j]] << std::endl;
    //     std::cout << "spider r theta " << angle->GetRecoilTheta(angle->GetProjectileSpiderTheta(ed->spider_ring_id[v_spider_p_index[j]], ed->spider_sector_id[v_spider_p_index[j]])) << std::endl;
    //     std::cout << "spider r phi " << angle->GetRecoilPhi(angle->GetProjectileSpiderPhi(ed->spider_ring_id[v_spider_p_index[j]], ed->spider_sector_id[v_spider_p_index[j]])) << std::endl;
    //   }
    // }
    // if(flag_s3_p){
    //   for(int j=0;j<num_s3_sector_p;j++){
    //     std::cout << "s3 sector p energy " << ed->s3_sector_energy[v_s3_sector_p_index[j]] << std::endl;
    //     std::cout << "s3 sector p phi " << angle->GetProjectileS3Phi(ed->s3_sector_id[v_s3_sector_p_index[j]]) << std::endl;
    //   }
    //   for(int j=0;j<num_s3_ring_p;j++){
    //     std::cout << "s3 ring p energy " << ed->s3_ring_energy[v_s3_ring_p_index[j]] << std::endl;
    //     std::cout << "s3 sector p theta " << angle->GetProjectileS3Theta(ed->s3_ring_id[v_s3_ring_p_index[j]]) << std::endl;
    //   }
    // }
    // if(flag_s3_r){
    //   for(int j=0;j<num_s3_sector_r;j++){
    //     std::cout << "s3 sector r energy " << ed->s3_sector_energy[v_s3_sector_r_index[j]] << std::endl;
    //     std::cout << "s3 sector r phi " << angle->GetRecoilPhi(angle->GetProjectileS3Theta(ed->s3_sector_id[v_s3_sector_p_index[j]])) << std::endl;
    //   }
    //   for(int j=0;j<num_s3_ring_p;j++){
    //     std::cout << "s3 ring r energy " << ed->s3_ring_energy[v_s3_ring_r_index[j]] << std::endl;
    //     std::cout << "s3 sector r theta " << angle->GetRecoilTheta(angle->GetProjectileS3Theta(ed->s3_ring_id[v_s3_ring_p_index[j]])) << std::endl;
    //   }
    // }
    // std::cout << "\n\n";

    // spider p
    if(flag_spider_p && !flag_spider_r && !flag_s3_p && !flag_s3_r){
      int ii = v_spider_p_index[0];
      
      for(size_t j=0;j<v_ge_index.size();j++){
        int jj = v_ge_index[j];
        rg[j] = ed->ge_ring_id[jj];
        sg[j] = ed->ge_sector_id[jj];

        e_raw[j] = ed->ge_energy[jj];

        tg = angle->GetGammaTheta(ed->ge_ring_id[jj], ed->ge_sector_id[jj]);
        pg = angle->GetGammaPhi(ed->ge_ring_id[jj], ed->ge_sector_id[jj]);

        tp_p = angle->GetProjectileSpiderTheta(ed->spider_ring_id[ii], ed->spider_sector_id[ii]);
        pp_p = angle->GetProjectileSpiderPhi(ed->spider_ring_id[ii], ed->spider_sector_id[ii]);
        cospg_p = angle->ComputeCosThetaParticleGamma(tg, pg, tp_p, pp_p);
        beta_p = angle->ComputeProjectileBeta(tp_p, 1000); // 1000fs means stopped
        e_dc_p[j] = ComputeEnergyWithDoppler(ed->ge_energy[jj], beta_p, cospg_p);

        tp_r = angle->GetRecoilTheta(tp_p);
        pp_r = angle->GetRecoilPhi(pp_p);
        cospg_r = angle->ComputeCosThetaParticleGamma(tg, pg, tp_r, pp_r);
        beta_r = angle->ComputeRecoilBeta(tp_r, 1000); // 1000fs means stopped
        e_dc_r[j] = ComputeEnergyWithDoppler(ed->ge_energy[jj], beta_r, cospg_r);

        // cout
        // std::cout << "ge ring " << rg[j] << " sector " << sg[j] << std::endl;   
        // std::cout << "ge ring " << ed->ge_ring_id[j] << " sector " << ed->ge_sector_id[j] << " energy " << ed->ge_energy[j] << " keV " << std::endl;
        // std::cout << "gamma angle tg " << tg << " pg " << pg << std::endl;
        // std::cout << "projectile tp_p " << tp_p << " pp_p " << pp_p << std::endl;
        // std::cout << "recoil tp_r " << tp_r << " pp_r " << pp_p << std::endl;
        // std::cout << "beta_p " << beta_p << " cospg_p " << cospg_p << std::endl;
        // std::cout << "beta_r " << beta_r << " cospg_r " << cospg_r << std::endl;
        // std::cout << "e_raw " << e_raw[j] << std::endl;
        // std::cout << "e_dc_p " << e_dc_p[j] << std::endl;
        // std::cout << "e_dc_r " << e_dc_r[j] << std::endl;
      }

      hits = v_ge_index.size();
      rp = ed->spider_ring_id[ii];
      sp = ed->spider_sector_id[ii];

      if(hits>0) tr->Fill();
      Clear();
    }

    // s3 p
    if(!flag_spider_p && !flag_spider_r && flag_s3_p && !flag_s3_r){
      int ii_sector = v_s3_sector_p_index[0];
      int ii_ring = v_s3_ring_p_index[0];

      for(size_t j=0;j<v_ge_index.size();j++){
        int jj = v_ge_index[j];
        rg[j] = ed->ge_ring_id[jj];
        sg[j] = ed->ge_sector_id[jj];

        e_raw[j] = ed->ge_energy[jj];

        tg = angle->GetGammaTheta(ed->ge_ring_id[jj], ed->ge_sector_id[jj]);
        pg = angle->GetGammaPhi(ed->ge_ring_id[jj], ed->ge_sector_id[jj]);

        tp_p = angle->GetProjectileS3Theta(ed->s3_ring_id[ii_ring]);
        pp_p = angle->GetProjectileS3Phi(ed->s3_sector_id[ii_sector]);
        cospg_p = angle->ComputeCosThetaParticleGamma(tg, pg, tp_p, pp_p);
        beta_p = angle->ComputeProjectileBeta(tp_p, 1000); // 1000fs means stopped
        e_dc_p[j] = ComputeEnergyWithDoppler(ed->ge_energy[jj], beta_p, cospg_p);

        tp_r = angle->GetRecoilTheta(tp_p);
        pp_r = angle->GetRecoilPhi(pp_p);
        cospg_r = angle->ComputeCosThetaParticleGamma(tg, pg, tp_r, pp_r);
        beta_r = angle->ComputeRecoilBeta(tp_r, 1000); // 1000fs means stopped
        e_dc_r[j] = ComputeEnergyWithDoppler(ed->ge_energy[jj], beta_r, cospg_r);

        // cout
        // std::cout << "ge ring " << rg[j] << " sector " << sg[j] << std::endl;   
        // std::cout << "ge ring " << ed->ge_ring_id[j] << " sector " << ed->ge_sector_id[j] << " energy " << ed->ge_energy[j] << " keV " << std::endl;
        // std::cout << "gamma tg " << tg << " pg " << pg << std::endl;
        // std::cout << "projectile tp_p " << tp_p << " pp_p " << pp_p << std::endl;
        // std::cout << "recoil tp_r " << tp_r << " pp_r " << pp_p << std::endl;
        // std::cout << "beta_p " << beta_p << " cospg_p " << cospg_p << std::endl;
        // std::cout << "beta_r " << beta_r << " cospg_r " << cospg_r << std::endl;
        // std::cout << "e_raw " << e_raw[j] << std::endl;
        // std::cout << "e_dc_p " << e_dc_p[j] << std::endl;
        // std::cout << "e_dc_r " << e_dc_r[j] << std::endl;
      }

      hits = v_ge_index.size();
      rp = ed->s3_ring_id[ii_ring];
      sp = ed->s3_sector_id[ii_sector];

      if(hits>0) tr->Fill();
      Clear();
    }

    // std::cout << "\nooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n";
  }
    
  std::cout << std::endl;
}

//
void doppler::SaveFile()
{
  file_out->cd();
  tr_event_out->Write();
  tr_bg_out->Write();
  file_out->Close();
}

//
double doppler::ComputeEnergyWithDoppler(double energy, double beta, double tg, double pg, double tp, double pp)
{
  tg = tg/180.*PI;
  pg = pg/180.*PI;
  tp = tp/180.*PI;
  pp = pp/180.*PI;

  double tmp1 = sin(tp)*sin(tg)*cos(pp-pg);
  double tmp2 = cos(tp)*cos(tg);

  double energy_doppler = energy*(1-beta*(tmp1+tmp2))/sqrt(1-beta*beta);

  return energy_doppler;
}

//
double doppler::ComputeEnergyWithDoppler(double energy, double beta, double cospg)
{
  double energy_doppler = energy*(1-beta*cospg)/sqrt(1-beta*beta);

  return energy_doppler; 
}

//
bool doppler::ReadCutInfo()
{
  std::cout << "read cut info" << std::endl;

  std::string line;
  int ring = 0, sector = 0;
  double min = 0, max = 0;
  int key = 0;

  // spider cut p
  std::ifstream fi_spider_cut_p;
  fi_spider_cut_p.open("../pars/cut_info/spider_cut_p.txt");
  if(!fi_spider_cut_p){
    std::cout << "can not open spider_cut_p.txt." << std::endl;
    return 0;
  }
  std::getline(fi_spider_cut_p, line);
  while(1){
    fi_spider_cut_p >> ring >> sector >> min >> max;
    if(!fi_spider_cut_p)  break;

    key = 100*ring + sector;
    map_spider_cut_p[key] = std::make_pair(min, max);
  }
  fi_spider_cut_p.close();

  // spider cut r
  std::ifstream fi_spider_cut_r;
  fi_spider_cut_r.open("../pars/cut_info/spider_cut_r.txt");
  if(!fi_spider_cut_r){
    std::cout << "can not open spider_cut_r.txt." << std::endl;
    return 0;
  }
  std::getline(fi_spider_cut_r, line);
  while(1){
    fi_spider_cut_r >> ring >> sector >> min >> max;
    if(!fi_spider_cut_r)  break;

    key = 100*ring + sector;
    map_spider_cut_r[key] = std::make_pair(min, max);
  }
  fi_spider_cut_r.close();

  // s3 sector cut p
  std::ifstream fi_s3_sector_cut_p;
  fi_s3_sector_cut_p.open("../pars/cut_info/s3_sector_cut_p.txt");
  if(!fi_s3_sector_cut_p){
    std::cout << "can not open s3_sector_cut_p.txt." << std::endl;
    return 0;
  }
  std::getline(fi_s3_sector_cut_p, line);
  while(1){
    fi_s3_sector_cut_p >> sector >> min >> max;
    if(!fi_s3_sector_cut_p)  break;

    key = sector;
    map_s3_sector_cut_p[key] = std::make_pair(min, max);
  }
  fi_s3_sector_cut_p.close();

  // s3 sector cut r
  std::ifstream fi_s3_sector_cut_r;
  fi_s3_sector_cut_r.open("../pars/cut_info/s3_sector_cut_r.txt");
  if(!fi_s3_sector_cut_r){
    std::cout << "can not open s3_sector_cut_r.txt." << std::endl;
    return 0;
  }
  std::getline(fi_s3_sector_cut_r, line);
  while(1){
    fi_s3_sector_cut_r >> sector >> min >> max;
    if(!fi_s3_sector_cut_r)  break;

    key = sector;
    map_s3_sector_cut_r[key] = std::make_pair(min, max);
  }
  fi_s3_sector_cut_p.close();

  // s3 ring cut p
  std::ifstream fi_s3_ring_cut_p;
  fi_s3_ring_cut_p.open("../pars/cut_info/s3_ring_cut_p.txt");
  if(!fi_s3_ring_cut_p){
    std::cout << "can not open s3_ring_cut_p.txt." << std::endl;
    return 0;
  }
  std::getline(fi_s3_ring_cut_p, line);
  while(1){
    fi_s3_ring_cut_p >> ring >> min >> max;
    if(!fi_s3_ring_cut_p)  break;

    key = ring;
    map_s3_ring_cut_p[key] = std::make_pair(min, max);
  }
  fi_s3_ring_cut_p.close();

  // s3 ring cut r
  std::ifstream fi_s3_ring_cut_r;
  fi_s3_ring_cut_r.open("../pars/cut_info/s3_ring_cut_r.txt");
  if(!fi_s3_ring_cut_r){
    std::cout << "can not open s3_ring_cut_r.txt." << std::endl;
    return 0;
  }
  std::getline(fi_s3_ring_cut_r, line);
  while(1){
    fi_s3_ring_cut_r >> ring >> min >> max;
    if(!fi_s3_ring_cut_r)  break;

    key = ring;
    map_s3_ring_cut_r[key] = std::make_pair(min, max);
  }
  fi_s3_ring_cut_r.close();

  return 1;
}

//
void doppler::PrintCutInfo()
{
  // spider cut p
  std::cout << "spider cut p cut info " << std::endl; 
  for(const auto &it : map_spider_cut_p){
    std::cout << "ring " << it.first/100 << " sector " << it.first%100 << " e_min " << it.second.first << " e_max " << it.second.second << std::endl;
  }

  // spider cut r
  std::cout << "spider cut r cut info " << std::endl; 
  for(const auto &it : map_spider_cut_r){
    std::cout << "ring " << it.first/100 << " sector " << it.first%100 << " e_min " << it.second.first << " e_max " << it.second.second << std::endl;
  }

  // s3 ring cut p
  std::cout << "s3 ring cut p cut info " << std::endl; 
  for(const auto &it : map_s3_ring_cut_p){
    std::cout << "ring " << it.first << " e_min " << it.second.first << " e_max " << it.second.second << std::endl;
  }

  // s3 ring cut r
  std::cout << "s3 ring cut r cut info " << std::endl; 
  for(const auto &it : map_s3_ring_cut_r){
    std::cout << "ring " << it.first << " e_min " << it.second.first << " e_max " << it.second.second << std::endl;
  }

  // s3 sector cut p
  std::cout << "s3 sector cut p cut info " << std::endl; 
  for(const auto &it : map_s3_sector_cut_p){
    std::cout << "sector " << it.first << " e_min " << it.second.first << " e_max " << it.second.second << std::endl;
  }

  // s3 sector cut r
  std::cout << "s3 sector cut r cut info " << std::endl; 
  for(const auto &it : map_s3_sector_cut_r){
    std::cout << "sector " << it.first << " e_min " << it.second.first << " e_max " << it.second.second << std::endl;
  }
}

//
bool doppler::ReadGeConfigInfo()
{
  std::cout << "read ge configure info" << std::endl;

  std::string line;
  int sid = 0, ch = 0;
  bool gb = 0;
  int key = 0;

  // ge configure
  std::ifstream fi;
  fi.open(TString::Format("../pars/config/%s.txt", configure.c_str()));
  if(!fi){
    std::cout << "can not open configure.txt." << std::endl;
    return 0;
  }
  std::getline(fi, line);
  while(1){
    fi >> sid >> ch >> gb;
    if(!fi)  break;

    key = 100*sid + ch;
    map_ge_flag[key] = gb;
  }
  fi.close();

  return 1;
}

//
void doppler::PrintGeConfigInfo()
{
  //
  std::cout << "ge configure info " << std::endl; 
  for(const auto &it : map_ge_flag){
    std::cout << "sid " << it.first/100 << " ch " << it.first%100 << " " << it.second << std::endl;
  }
}

//
void doppler::Clear()
{
  hits = 0;
  flag_spider_p = 0; flag_spider_r = 0;
  flag_s3_p = 0; flag_s3_r = 0;
  rp = 0; sp = 0;
  rr = 0; sr = 0;
  memset(rg, 0, sizeof(rg));
  memset(sg, 0, sizeof(sg));

  memset(e_raw, 0, sizeof(e_raw));
  memset(e_dc_p, 0, sizeof(e_dc_p));
  memset(e_dc_r, 0, sizeof(e_dc_r));
}

// Current build processing and parameter scans share Model::evaluate above.
namespace sm150_doppler {
namespace {
using json=nlohmann::json;
namespace fs=std::filesystem;
void require(bool ok,const std::string &message) { if(!ok)throw std::runtime_error(message); }
json read_json(const std::string &path) {
  std::ifstream in(path);require(bool(in),"Cannot open "+path);json j;in>>j;return j;
}
json object(TFile &f,const char *key) {
  auto *x=dynamic_cast<TObjString*>(f.Get(key));require(x,std::string("Missing object ")+key);
  return json::parse(x->GetString().Data());
}
std::string sha256(const fs::path &path) {
  std::ifstream in(path,std::ios::binary);require(bool(in),"Cannot hash "+path.string());
  std::unique_ptr<EVP_MD_CTX,decltype(&EVP_MD_CTX_free)> ctx(EVP_MD_CTX_new(),EVP_MD_CTX_free);
  require(ctx&&EVP_DigestInit_ex(ctx.get(),EVP_sha256(),nullptr)==1,"SHA initialization failed");
  std::vector<char> data(1024*1024);
  while(in) { in.read(data.data(),data.size());require(EVP_DigestUpdate(ctx.get(),data.data(),in.gcount())==1,"SHA update failed"); }
  require(in.eof(),"SHA read failed");unsigned char bytes[EVP_MAX_MD_SIZE];unsigned int size;
  require(EVP_DigestFinal_ex(ctx.get(),bytes,&size)==1,"SHA finish failed");
  std::ostringstream out;for(unsigned int i=0;i<size;++i)out<<std::hex<<std::setw(2)<<std::setfill('0')<<int(bytes[i]);return out.str();
}
std::string period(int group) { return std::array<std::string,3>{"high1","low","high2"}.at(group); }
int group_id(const std::string &group) {
  for(int i=0;i<3;++i)if(period(i)==group)return i;
  throw std::runtime_error("Unknown build period");
}
struct Histograms {
  std::map<std::string,std::unique_ptr<TH1D>> data;
  int bins;double low,high;
  Histograms(int n,double a,double b):bins(n),low(a),high(b){}
  void fill(const std::string &name,double energy,double weight) {
    auto &h=data[name];if(!h) {
      h=std::make_unique<TH1D>(name.c_str(),";Gamma energy (keV);Weighted pairs",bins,low,high);
      h->SetDirectory(nullptr);h->Sumw2();
    }
    h->Fill(energy,weight);
  }
  void write(TFile &f) {f.cd();for(auto &h:data)require(h.second->Write()>0,"Histogram write failed");}
};
struct ParticleReader {
  TTreeReaderArray<Short_t> sid,ch,side,id;
  std::unique_ptr<TTreeReaderArray<Short_t>> sector;
  TTreeReaderArray<Double_t> energy,dt,alpha;
  TTreeReaderArray<Long64_t> source;
  ParticleReader(TTreeReader &r,const std::string &name,const std::string &index):
    sid(r,(name+"_sid").c_str()),ch(r,(name+"_ch").c_str()),side(r,(name+"_random_side").c_str()),
    id(r,(name+"_"+index).c_str()),energy(r,(name+"_energy").c_str()),dt(r,(name+"_dt_ns").c_str()),
    alpha(r,(name+"_alpha").c_str()),source(r,(name+"_source_entry").c_str()) {
    if(name=="spider")sector=std::make_unique<TTreeReaderArray<Short_t>>(r,"spider_sector_id");
  }
  size_t size()const{return energy.GetSize();}
  void validate()const {
    const auto n=size();require(sid.GetSize()==n&&ch.GetSize()==n&&side.GetSize()==n&&id.GetSize()==n&&dt.GetSize()==n&&alpha.GetSize()==n&&source.GetSize()==n&&(!sector||sector->GetSize()==n),"Inconsistent particle arrays");
  }
};
struct PairReader {
  TTreeReader r;
  TTreeReaderValue<Int_t> run,group,ge,gr,gs,ref,particle,ring,sector,side,flags;
  TTreeReaderValue<Long64_t> ge_source,particle_source,sector_source;
  TTreeReaderValue<Double_t> raw,energy,dt,alpha,extra;
  TTreeReaderValue<Bool_t> legacy;
  explicit PairReader(TTree *t):r(t),run(r,"run"),group(r,"group"),ge(r,"ge"),gr(r,"ge_ring"),gs(r,"ge_sector"),
    ref(r,"reference"),particle(r,"particle"),ring(r,"ring"),sector(r,"sector"),side(r,"random_side"),flags(r,"ge_flags"),
    ge_source(r,"ge_source_entry"),particle_source(r,"particle_source_entry"),sector_source(r,"sector_source_entry"),
    raw(r,"e_raw"),energy(r,"particle_MeV"),dt(r,"dt_ns"),alpha(r,"alpha"),extra(r,"extra_factor"),legacy(r,"exclusive_max_energy"){}
  Pair get() {
    Pair p;p.run=*run;p.group=*group;p.ge=*ge;p.ge_ring=*gr;p.ge_sector=*gs;p.reference=*ref;
    p.particle=*particle;p.ring=*ring;p.sector=*sector;p.side=*side;p.ge_flags=*flags;
    p.ge_source=*ge_source;p.particle_source=*particle_source;p.sector_source=*sector_source;
    p.raw=*raw;p.particle_MeV=*energy;p.dt_ns=*dt;p.alpha=*alpha;p.extra_factor=*extra;p.legacy_selection=*legacy;return p;
  }
};
json provenance(const Options &o) {
  json hashes=json::object();
  for(const auto &path:{"angle_info/ge_angle.txt","angle_info/spider_angle.txt","angle_info/s3_ring_angle.txt","angle_info/s3_sector_angle.txt","beta_info/theta_ein_eout_p.txt","beta_info/theta_ein_eout_r.txt"})
    hashes[path]=sha256(fs::path(o.project)/"pars"/path);
  return {{"angle_and_beta_sha256",hashes},{"executable_sha256",sha256("/proc/self/exe")},
    {"ROOT_version",gROOT->GetVersion()},{"mass_convention","Historical mass expressions retained; not new nuclear-mass calibration"}};
}
void finish(TFile &out,const json &manifest) {
  out.cd();TObjString s(manifest.dump().c_str());require(s.Write("sm150_doppler_manifest")>0,"Manifest write failed");out.Close();
}
}
int process_build(const Options &o) {
  require(!o.input.empty()&&!o.parameters.empty(),"Input and parameters are required");
  require(!fs::exists(o.output),"Output exists; refusing overwrite");
  auto parameters=read_json(o.parameters);Model model(o.project,parameters);
  const auto source_sha=sha256(o.input);
  const auto input_size=fs::file_size(o.input);
  const auto input_time=fs::last_write_time(o.input);
  std::unique_ptr<TFile> in(TFile::Open(o.input.c_str(),"READ"));require(in&&!in->IsZombie(),"Invalid build file");
  const auto m=object(*in,"sm150_build_manifest"),expected=object(*in,"sm150_build_counts");
  require(m.at("schema_version")==2&&m.at("selection").at("single_Ge")==true,"Require schema-2 single-Ge build");
  const int run=std::stoi(m.at("run").get<std::string>()),group=group_id(m.at("group"));
  const auto receipt=read_json((fs::path(o.input).parent_path()/("data"+m.at("run").get<std::string>()+"_build_cpp_complete.json")).string());
  bool hash_ok=false;for(const auto &x:receipt.at("outputs"))if(fs::path(x.at("path").get<std::string>()).filename()==fs::path(o.input).filename())hash_ok=x.at("sha256")==source_sha&&x.at("bytes")==input_size;
  require(hash_ok,"Build hash or size disagrees with completion receipt");
  const fs::path policy_path=fs::path(o.project)/"doppler/parameters/random_subtraction.json";
  const auto policy=read_json(policy_path.string());
  std::map<std::pair<int,int>,json> gates;for(const auto &g:m.at("windows"))gates.emplace(std::make_pair(g.at("ge").get<int>(),g.at("reference").get<int>()),g);
  TFile out(o.output.c_str(),"NEW","Pair-preserving Doppler correction",404);require(!out.IsZombie(),"Cannot create output");
  TTree tree("pairs","Accepted build pairs; no additional timing or energy cuts");
  Pair p;Correction c;double weight;
  auto integer=[&](const char *n,int *v){tree.Branch(n,v,(std::string(n)+"/I").c_str());};
  integer("run",&p.run);integer("group",&p.group);integer("ge",&p.ge);integer("ge_ring",&p.ge_ring);integer("ge_sector",&p.ge_sector);
  integer("reference",&p.reference);integer("particle",&p.particle);integer("ring",&p.ring);integer("sector",&p.sector);integer("random_side",&p.side);integer("ge_flags",&p.ge_flags);
  tree.Branch("ge_source_entry",&p.ge_source,"ge_source_entry/L");tree.Branch("particle_source_entry",&p.particle_source,"particle_source_entry/L");tree.Branch("sector_source_entry",&p.sector_source,"sector_source_entry/L");
  auto real=[&](const char *n,double *v){tree.Branch(n,v,(std::string(n)+"/D").c_str());};
  real("e_raw",&p.raw);real("particle_MeV",&p.particle_MeV);real("dt_ns",&p.dt_ns);real("alpha",&p.alpha);real("extra_factor",&p.extra_factor);real("weight",&weight);
  real("e_dc_r",&c.energy);real("beta",&c.beta);real("cos_psi",&c.cosine);real("theta_p_deg",&c.theta_p);real("phi_p_deg",&c.phi_p);real("theta_r_deg",&c.theta_r);real("phi_r_deg",&c.phi_r);real("theta_ge_deg",&c.theta_g);real("phi_ge_deg",&c.phi_g);real("expected_particle_MeV",&c.expected_particle_MeV);
  tree.Branch("exclusive_max_energy",&p.legacy_selection,"exclusive_max_energy/O");tree.Branch("kinematic_extrapolated",&c.extrapolated,"kinematic_extrapolated/O");tree.Branch("stopped_candidate",&c.stopped,"stopped_candidate/O");
  Histograms hist(2200,0,2200);json counts=json::object();
  for(const auto &name:{"tr_event","tr_bg"}) {
    auto *t=dynamic_cast<TTree*>(in->Get(name));require(t,"Missing build tree");const bool prompt=std::string(name)=="tr_event";
    TTreeReader reader(t);TTreeReaderValue<Int_t> input_run(reader,"run"),n_ge(reader,"n_ge");
    TTreeReaderArray<Short_t> sid(reader,"ge_sid"),ch(reader,"ge_ch"),gr(reader,"ge_ring_id"),gs(reader,"ge_sector_id");
    TTreeReaderArray<Double_t> e(reader,"ge_energy");TTreeReaderArray<UChar_t> flags(reader,"ge_flags");TTreeReaderArray<Long64_t> source(reader,"ge_source_entry");
    ParticleReader si(reader,"spider","ring_id"),ring(reader,"s3_ring","id"),sector(reader,"s3_sector","id");
    Long64_t rows=0,ns=0,nr=0,legacy=0;
    while(reader.Next()) {
      ++rows;require(*n_ge==1&&*input_run==run&&sid.GetSize()==1&&ch.GetSize()==1&&gr.GetSize()==1&&gs.GetSize()==1&&e.GetSize()==1&&flags.GetSize()==1&&source.GetSize()==1,"Invalid Ge row");
      si.validate();ring.validate();sector.validate();require(ring.size()==sector.size()&&ring.size()<=1,"Invalid S3 pairing");
      if(ring.size()) require(ring.side[0]==sector.side[0]&&ring.alpha[0]==sector.alpha[0],"S3 weight/side mismatch");
      size_t max_si=0;for(size_t i=1;i<si.size();++i)if(si.energy[i]>si.energy[max_si])max_si=i;
      for(int ref=0;ref<2;++ref) {
        const auto &part=ref ? ring:si;
        for(size_t i=0;i<part.size();++i) {
          p=Pair{};p.run=run;p.group=group;p.ge=100*sid[0]+ch[0];p.ge_ring=gr[0];p.ge_sector=gs[0];p.reference=ref;
          p.particle=10000+100*part.sid[i]+part.ch[i];p.ring=part.id[i];p.sector=ref ? sector.id[0]:(*si.sector)[i];
          p.side=part.side[i];p.ge_flags=flags[0];p.ge_source=source[0];p.particle_source=part.source[i];p.sector_source=ref ? sector.source[0]:-1;
          p.raw=e[0];p.particle_MeV=part.energy[i]/1000;p.dt_ns=part.dt[i];p.alpha=part.alpha[i];
          p.extra_factor=policy.at("alpha_by_group").at((group==1 ? std::string("low_"):std::string("high_"))+(ref ? "s3":"si"));
          require(prompt ? p.side==0 : (p.side==-1||p.side==1),"Wrong prompt/random side");
          const auto &gate=gates.at({p.ge,ref});require(std::abs(p.alpha-gate.at("alpha").get<double>())<1e-12,"Stored pair alpha mismatch");
          bool accepted=prompt ? p.dt_ns>=gate.at("lower_ns").get<double>()&&p.dt_ns<gate.at("upper_ns").get<double>() : false;
          if(!prompt)for(const auto &range:gate.at("random"))if(p.dt_ns>=range[0].get<double>()&&p.dt_ns<range[1].get<double>())accepted=true;
          require(accepted,"Pair outside its embedded gate");
          p.legacy_selection=ref ? si.size()==0 : (ring.size()==0&&i==max_si);
          c=model.evaluate(p);weight=prompt ? 1 : -p.alpha*p.extra_factor;
          tree.Fill();if(ref)++nr;else ++ns;if(p.legacy_selection)++legacy;
          const std::string family=ref ? "S3":"Si";
          for(const auto &spec:std::array<std::pair<std::string,double>,2>{{{"raw",p.raw},{"corrected",c.energy}}}) {
            const auto base=family+"_"+spec.first;
            hist.fill(base+"_net",spec.second,weight);
            hist.fill(base+(prompt ? "_prompt":"_random_scaled"),spec.second,prompt ? 1 : -weight);
            if(!prompt)hist.fill(base+"_random_raw",spec.second,1);
            if(p.legacy_selection)hist.fill(base+"_exclusive_max_energy_net",spec.second,weight);
          }
        }
      }
    }
    require(reader.GetEntryStatus()==TTreeReader::kEntryBeyondEnd,"Tree reading failed");
    require(rows==expected.at(name).at("events")&&ns==expected.at(name).at("si_pairs")&&nr==expected.at(name).at("s3_pairs"),"Input pair counts do not match build metadata");
    counts[name]={{"events",rows},{"si_pairs",ns},{"s3_pairs",nr},{"exclusive_max_energy_pairs",legacy}};
  }
  require(input_size==fs::file_size(o.input)&&input_time==fs::last_write_time(o.input),"Input changed during processing");
  out.cd();require(tree.Write()>0,"Pair tree write failed");hist.write(out);
  json manifest=provenance(o);
  manifest.update({{"schema_version",1},{"mode","current_build_pairs"},{"source",fs::absolute(o.input).string()},{"source_sha256",source_sha},{"build_manifest",m},{"parameters",parameters},{"parameters_sha256",sha256(o.parameters)},{"counts",counts},{"pair_entries",tree.GetEntries()},{"random_policy",policy},{"random_policy_sha256",sha256(policy_path)},{"exclusive_max_energy_definition","Exclusive Si/S3 record; highest-energy Si only. A selection diagnostic, not a complete legacy replay."}});
  finish(out,manifest);std::cout<<json({{"output",o.output},{"counts",counts}}).dump()<<"\n";return 0;
}
int scan_pairs(const Options &o) {
  require(!o.scan_parameters.empty()&&!fs::exists(o.output),"Scan parameters missing or output exists");
  const auto configs=read_json(o.scan_parameters);require(configs.is_array()&&!configs.empty()&&configs.size()<=100,"Expect 1-100 scan configurations");
  std::vector<Model> models;for(const auto &c:configs)models.emplace_back(o.project,c);
  Histograms hist(1040,280,800);json inputs=json::array();std::vector<Long64_t> extrapolated(configs.size()),stopped(configs.size());
  Long64_t pairs=0;std::ifstream list(o.scan_list);require(bool(list),"Cannot read scan input list");std::string path;
  while(std::getline(list,path)) {
    if(path.empty())continue;
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(),"READ"));require(f&&!f->IsZombie(),"Invalid pair file");
    const auto manifest=object(*f,"sm150_doppler_manifest");require(manifest.at("mode")=="current_build_pairs","Not a pair-preserving Doppler file");
    auto *tree=dynamic_cast<TTree*>(f->Get("pairs"));require(tree,"Missing pair tree");
    inputs.push_back({{"path",path},{"sha256",sha256(path)},{"pair_entries",tree->GetEntries()},{"build_source_sha256",manifest.at("source_sha256")}});
    PairReader reader(tree);Long64_t n=0;
    while(reader.r.Next()) {
      const auto p=reader.get();++pairs;++n;const double w=p.side ? -p.alpha*p.extra_factor : 1;
      const std::string suffix=period(p.group)+(p.reference ? "_S3":"_Si");
      for(size_t i=0;i<models.size();++i) {
        auto c=models[i].evaluate(p);extrapolated[i]+=c.extrapolated;stopped[i]+=c.stopped;
        const std::string prefix="v"+std::to_string(i)+"_"+suffix;
        hist.fill(prefix+"_all",c.energy,w);hist.fill(prefix+"_ge"+std::to_string(p.ge),c.energy,w);
      }
    }
    require(reader.r.GetEntryStatus()==TTreeReader::kEntryBeyondEnd&&n==tree->GetEntries(),"Incomplete pair scan");
    std::cout<<"Scanned "<<path<<" pairs="<<n<<std::endl;
  }
  require(pairs>0,"No scan input pairs");TFile out(o.output.c_str(),"NEW","Doppler scan",404);require(!out.IsZombie(),"Cannot create scan output");hist.write(out);
  auto manifest=provenance(o);manifest.update({{"schema_version",1},{"mode","parameter_scan"},{"parameters",configs},{"input_files",inputs},{"pairs",pairs},{"extrapolated",extrapolated},{"stopped",stopped},{"histogram_bin_keV",0.5},{"histogram_range_keV",{280,800}},{"variance","sum of squared signed pair weights"}});
  finish(out,manifest);std::cout<<json({{"output",o.output},{"pairs",pairs},{"variants",configs.size()}}).dump()<<"\n";return 0;
}
}
