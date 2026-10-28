#include "build.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include "TNamed.h"

//
build::build(const std::string &filename_in, const std::string &filename_out, int r)
{
  benchmark = new TBenchmark;

  run = r;
  if((run>=457 && run<=462) || (run>=472 && run<=605) || (run>=680 && run<=718)){
    v_s3_sid = {9, 10, 11, 12};
  }else if((run>=616 && run<=674) || (run>=721 && run<=723)){
    v_s3_sid = {2, 3, 11, 12};
  }else{
    throw std::invalid_argument("wrong run number.");
  }

  if(!InitMapSectorRingID()){
    throw std::invalid_argument("can not init ring and sector id data.");
  }
  PrintMapSectorRingID();

  file_in = TFile::Open(filename_in.c_str());
  if(file_in->IsZombie()){
    std::cout << "open file " << filename_in << " error!" << std::endl;
    delete file_in;

    return ;
  }

  tr_in = (TTree*)file_in->Get("tree");
  tr_in->SetBranchAddress("sr", &sr);
  tr_in->SetBranchAddress("cid", &cid);
  tr_in->SetBranchAddress("sid", &sid);
  tr_in->SetBranchAddress("ch", &ch);
  tr_in->SetBranchAddress("evte", &evte);
  tr_in->SetBranchAddress("ts", &ts);
  tr_in->SetBranchAddress("cfd", &cfd);
  tr_in->SetBranchAddress("cfdft", &cfdft);
  tr_in->SetBranchAddress("cfds", &cfds);


  //
  if(std::strcmp(TF,"tf1") == 0) tf = new TF1("tf", "[0]+[1]/pow(x,0.5)+[2]/pow(x,1)+[3]/pow(x,2)", CUTGE, 4096);
  else if(std::strcmp(TF, "tf2") == 0) tf = new TF1("tf","[0]+[1]*exp([2]/pow(x,1.5))", CUTGE, 4096);
  else{
    throw std::invalid_argument("wrong tf name.");
  }

  // read cali data
  if(!ReadCaliData()){
    throw std::invalid_argument("can not read cali data.");
  }
  PrintCaliData();

  // read s3 cor
  ReadS3CorData();
  PrintS3CorData();

  time_alignment.Load("../pars", run);
  if(!time_alignment.total()){
    if(!ReadTsOffset()){
      throw std::invalid_argument("can not read ts offset data.");
    }
    PrintTsOffset();
    time_alignment.ValidateNativeOffsets(map_ts_offset);
  }
  std::cout << "time alignment profile: " << time_alignment.profile() << std::endl;

  // read ge time walk
  if(!time_alignment.active()){
    if(!ReadGeTwData()){
      throw std::invalid_argument("can not read ge time walk data.");
    }
    PrintGeTwData();
  }

  rndm = new TRandom3((Long64_t)time(0));
  file_out = TFile::Open(filename_out.c_str(), "recreate");
  file_out->cd();
  TNamed("sm150_time_alignment_profile", time_alignment.profile().c_str()).Write();
  TNamed("sm150_time_alignment_source_sha256", time_alignment.source_sha256().c_str()).Write();
}

//
build::~build()
{
  delete file_in;
  delete file_out;
}

//
void build::Process()
{
  benchmark->Start("build");

  //
  GetGeSpiderS3EventPrompt();
  GetGeSpiderS3EventRandom();
  
  benchmark->Show("build");
}

//
void build::GetGeSpiderS3EventPrompt()
{
  std::cout << "start get ge si prompt events" << std::endl;

  int n_max_ge = GENUM;
  int n_max_spider = SPIDERNUM;
  int n_max_s3_sector = S3SECTORNUM;
  int n_max_s3_ring = S3RINGNUM;

  double cut_ge = CUTGE; 
  double cut_si = CUTSI;

  double t_win = TIMEWINDOW;

  //
  Int_t n_ge = 0;
  Short_t ge_sid[n_max_ge];
  Short_t ge_ch[n_max_ge];
  Double_t ge_energy[n_max_ge];
  Long64_t ge_ts[n_max_ge];

  Int_t n_spider = 0;
  Short_t spider_sid[n_max_spider];
  Short_t spider_ch[n_max_spider];
  Double_t spider_energy[n_max_spider];
  Long64_t spider_ts[n_max_spider];

  Int_t n_s3_sector = 0;
  Short_t s3_sector_sid[n_max_s3_sector];
  Short_t s3_sector_ch[n_max_s3_sector];
  Short_t s3_sector_id[n_max_s3_sector];
  Double_t s3_sector_energy[n_max_s3_sector];
  Long64_t s3_sector_ts[n_max_s3_sector];
  Int_t n_s3_ring = 0;
  Short_t s3_ring_sid[n_max_s3_ring];
  Short_t s3_ring_ch[n_max_s3_ring];
  Short_t s3_ring_id[n_max_s3_ring];
  Double_t s3_ring_energy[n_max_s3_ring];
  Long64_t s3_ring_ts[n_max_s3_ring];

  std::stringstream ss;
  ss.str("");
  ss << "coincidence window [-" << t_win << ", " << t_win << "] ns";
  TTree *tr = new TTree("tr_event", ss.str().c_str());

  //
  memset(ge_sid, 0, sizeof(ge_sid));
  memset(ge_ch, 0, sizeof(ge_ch));
  memset(ge_energy, 0, sizeof(ge_energy));
  memset(ge_ts, 0, sizeof(ge_ts));

  memset(spider_sid, 0, sizeof(spider_sid));
  memset(spider_ch, 0, sizeof(spider_ch));
  memset(spider_energy, 0, sizeof(spider_energy));
  memset(spider_ts, 0, sizeof(spider_ts));

  memset(s3_sector_sid, 0, sizeof(s3_sector_sid));
  memset(s3_sector_ch, 0, sizeof(s3_sector_ch));
  memset(s3_sector_id, 0, sizeof(s3_sector_id));
  memset(s3_sector_energy, 0, sizeof(s3_sector_energy));
  memset(s3_sector_ts, 0, sizeof(s3_sector_ts));

  memset(s3_ring_sid, 0, sizeof(s3_ring_sid));
  memset(s3_ring_ch, 0, sizeof(s3_ring_ch));
  memset(s3_ring_id, 0, sizeof(s3_ring_id));
  memset(s3_ring_energy, 0, sizeof(s3_ring_energy));
  memset(s3_ring_ts, 0, sizeof(s3_ring_ts));

  //
  tr->Branch("n_ge", &n_ge, "n_ge/I");
  tr->Branch("ge_sid", ge_sid, "ge_sid[n_ge]/S");
  tr->Branch("ge_ch", ge_ch, "ge_ch[n_ge]/S");
  tr->Branch("ge_energy", ge_energy, "ge_energy[n_ge]/D");
  tr->Branch("ge_ts", ge_ts, "ge_ts[n_ge]/L");

  tr->Branch("n_spider", &n_spider, "n_spider/I");
  tr->Branch("spider_sid", spider_sid, "spider_sid[n_spider]/S");
  tr->Branch("spider_ch", spider_ch, "spider_ch[n_spider]/S");
  tr->Branch("spider_energy", spider_energy, "spider_energy[n_spider]/D");
  tr->Branch("spider_ts", spider_ts, "spider_ts[n_spider]/L");
  
  tr->Branch("n_s3_sector", &n_s3_sector, "n_s3_sector/I");
  tr->Branch("s3_sector_sid", s3_sector_sid, "s3_sector_sid[n_s3_sector]/S");
  tr->Branch("s3_sector_ch", s3_sector_ch, "s3_sector_ch[n_s3_sector]/S");
  tr->Branch("s3_sector_id", s3_sector_id, "s3_sector_id[n_s3_sector]/S");
  tr->Branch("s3_sector_energy", s3_sector_energy, "s3_sector_energy[n_s3_sector]/D");
  tr->Branch("s3_sector_ts", s3_sector_ts, "s3_sector_ts[n_s3_sector]/L");

  tr->Branch("n_s3_ring", &n_s3_ring, "n_s3_ring/I");
  tr->Branch("s3_ring_sid", s3_ring_sid, "s3_ring_sid[n_s3_ring]/S");
  tr->Branch("s3_ring_ch", s3_ring_ch, "s3_ring_ch[n_s3_ring]/S");
  tr->Branch("s3_ring_id", s3_ring_id, "s3_ring_id[n_s3_ring]/S");
  tr->Branch("s3_ring_energy", s3_ring_energy, "s3_ring_energy[n_s3_ring]/D");
  tr->Branch("s3_ring_ts", s3_ring_ts, "s3_ring_ts[n_s3_ring]/L");
  
  //
  Long64_t i_start = 0;

  Short_t sid1 = 0;
  Short_t ch1 = 0;
  Double_t energy1 = 0;
  Long64_t ts1 = 0;
  
  Short_t sid2 = 0;
  Short_t ch2 = 0;
  Double_t energy2 = 0;
  Long64_t ts2 = 0;

  while(true){//get first Ge data
    tr_in->GetEntry(i_start);
    energy = GetEnergy();
    ts_ns = GetTsns();
    i_start++;
    if(cid==0 && energy>cut_ge){
      break;
    }
  }
  
  Long64_t n_evt = 0;
  Long64_t i = i_start;
  Long64_t i_current = i-1;
  while(true){
    if(i==tr_in->GetEntries()) break;

    tr_in->GetEntry(i_current);
    energy = GetEnergy();
    ts_ns = GetTsns();

    sid1 = sid;
    ch1 = ch;
    energy1 = energy;
    ts1 = ts_ns;

    ge_sid[n_ge] = sid1;
    ge_ch[n_ge] = ch1; 
    ge_energy[n_ge] = energy1; 
    ge_ts[n_ge] = ts1; 
    n_ge++;

    // std::cout << std::endl;
    // std::cout << "info ge " << cid << " " << sid << " " << ch << " " << energy << " " << ts << std::endl;

    while(true){//search backward
      if(i_current--==0) break;
      tr_in->GetEntry(i_current);
      energy = GetEnergy();
      ts_ns = GetTsns();
      // std::cout << "info backward " << cid << " " << sid << " " << ch << " " << energy << " " << ts << std::endl;
      // std::cout << "time diff " << abs(ts-ts1) << std::endl;

      if(cid==0 && energy<=cut_ge){//if Ge data and small energy
        continue;
      }
      if(cid==1 && energy<=cut_si){//if Si/S3 data and small energy
        continue;
      }

      sid2 = sid;
      ch2 = ch;
      energy2 = energy;
      ts2 = ts_ns;
    
      if((abs(ts2-ts1))<=t_win){//this coincidence
        if(cid==0){
          ge_sid[n_ge] = sid2; 
          ge_ch[n_ge] = ch2; 
          ge_energy[n_ge] = energy2; 
          ge_ts[n_ge] = ts2; 
          n_ge++;
        }
        if(cid==1){
          if(!(std::find(v_s3_sid.begin(),v_s3_sid.end(),sid) != v_s3_sid.end())){// spider data
            spider_sid[n_spider] = sid2;
            spider_ch[n_spider] = ch2;
            spider_energy[n_spider] = energy2;
            spider_ts[n_spider] = ts2;
            n_spider++;
          }else{
            if(sid>=11&&sid<=12){
              s3_ring_sid[n_s3_ring] = sid2;
              s3_ring_ch[n_s3_ring] = ch2;
              s3_ring_id[n_s3_ring] = map_s3_ring_id[10000*cid+100*sid+ch];
              s3_ring_energy[n_s3_ring] = energy2;
              s3_ring_ts[n_s3_ring] = ts2;
              n_s3_ring++;
            }else{
              s3_sector_sid[n_s3_sector] = sid2;
              s3_sector_ch[n_s3_sector] = ch2;
              s3_sector_id[n_s3_sector] = map_s3_sector_id[10000*cid+100*sid+ch];
              s3_sector_energy[n_s3_sector] = energy2;
              s3_sector_ts[n_s3_sector] = ts2;
              n_s3_sector++;
            }
          }
        }
      }else break;
    }//while search backward

    while(true){//search forward
      if(i>=tr_in->GetEntries()) break;
      tr_in->GetEntry(i);
      energy = GetEnergy();
      ts_ns = GetTsns();
      // std::cout << "info forward " << cid << " " << sid << " " << ch << " " << energy << " " << ts << std::endl;
      // std::cout << "time diff " << abs(ts-ts1) << std::endl;

      if(cid==0 && energy<=cut_ge){//if ge data and small energy
        i++;
        continue;
      }
      if(cid==1 && energy<=cut_si){//if si data and small energy
        i++;
        continue;
      }

      sid2 = sid;
      ch2 = ch;
      energy2 = energy;
      ts2 = ts_ns;
    
      if((abs(ts2-ts1))<=t_win){//this coincidence
        i++;
        if(cid==0){
          ge_sid[n_ge] = sid2; 
          ge_ch[n_ge] = ch2; 
          ge_energy[n_ge] = energy2; 
          ge_ts[n_ge] = ts2; 
          n_ge++;
        }
        if(cid==1){
          if(!(std::find(v_s3_sid.begin(),v_s3_sid.end(),sid) != v_s3_sid.end())){// spider data
            spider_sid[n_spider] = sid2;
            spider_ch[n_spider] = ch2;
            spider_energy[n_spider] = energy2;
            spider_ts[n_spider] = ts2;
            n_spider++;
          }else{
            if(sid>=11&&sid<=12){
              s3_ring_sid[n_s3_ring] = sid2;
              s3_ring_ch[n_s3_ring] = ch2;
              s3_ring_id[n_s3_ring] = map_s3_ring_id[10000*cid+100*sid+ch];
              s3_ring_energy[n_s3_ring] = energy2;
              s3_ring_ts[n_s3_ring] = ts2;
              n_s3_ring++;
            }else{
              s3_sector_sid[n_s3_sector] = sid2;
              s3_sector_ch[n_s3_sector] = ch2;
              s3_sector_id[n_s3_sector] = map_s3_sector_id[10000*cid+100*sid+ch];
              s3_sector_energy[n_s3_sector] = energy2;
              s3_sector_ts[n_s3_sector] = ts2;
              n_s3_sector++;
            }
          }
        }
      }else break;
    } // while search forward

    // std::cout << "n_ge " << n_ge << std::endl;
    // std::cout << " n_spider " << n_spider << std::endl;
    // std::cout << " n_s3_sector " << n_s3_sector << " n_s3_ring " << n_s3_ring << std::endl;
    if(n_ge*(n_spider+n_s3_sector+n_s3_ring) > 0){
      n_evt++;

      file_out->cd();
      tr->Fill();

      if(n_evt%10000==0){
        std::cout << "\r" << n_evt << "  " << i << "/" << tr_in->GetEntries();
        std::cout << std::flush;
      }
    }

    n_ge = 0;
    memset(ge_sid, 0, sizeof(ge_sid));
    memset(ge_ch, 0, sizeof(ge_ch));
    memset(ge_energy, 0, sizeof(ge_energy));
    memset(ge_ts, 0, sizeof(ge_ts));

    n_spider = 0;
    memset(spider_sid, 0, sizeof(spider_sid));
    memset(spider_ch, 0, sizeof(spider_ch));
    memset(spider_energy, 0, sizeof(spider_energy));
    memset(spider_ts, 0, sizeof(spider_ts));

    n_s3_sector = 0;
    memset(s3_sector_sid, 0, sizeof(s3_sector_sid));
    memset(s3_sector_ch, 0, sizeof(s3_sector_ch));
    memset(s3_sector_id, 0, sizeof(s3_sector_id));
    memset(s3_sector_energy, 0, sizeof(s3_sector_energy));
    memset(s3_sector_ts, 0, sizeof(s3_sector_ts));

    n_s3_ring = 0;
    memset(s3_ring_sid, 0, sizeof(s3_ring_sid));
    memset(s3_ring_ch, 0, sizeof(s3_ring_ch));
    memset(s3_ring_id, 0, sizeof(s3_ring_id));
    memset(s3_ring_energy, 0, sizeof(s3_ring_energy));
    memset(s3_ring_ts, 0, sizeof(s3_ring_ts));

    while(true){//get next Ge data
      if(i>=tr_in->GetEntries()) break;
      tr_in->GetEntry(i);
      energy = GetEnergy();
      ts_ns = GetTsns();
      i++;
      if(cid==0 && energy>cut_ge){
        break;
      }
    }
    i_current = i-1;

  }//while

  std::cout << std::endl;

  file_out->cd();
  tr->Write();
}


//
void build::GetGeSpiderS3EventRandom()
{
  std::cout << "start get ge si random events" << std::endl;

  int n_max_ge = GENUM;
  int n_max_spider = SPIDERNUM;
  int n_max_s3_sector = S3SECTORNUM;
  int n_max_s3_ring = S3RINGNUM;

  double cut_ge = CUTGE; 
  double cut_si = CUTSI;

  double t_prompt = (TIMEWINDOW>TIMEJUMP) ? TIMEJUMP : TIMEWINDOW;
  double t_min_rand = TIMEJUMP;
  double t_max_rand = TIMEJUMP+TIMEWINDOW;

  //
  Int_t n_ge = 0;
  Short_t ge_sid[n_max_ge];
  Short_t ge_ch[n_max_ge];
  Double_t ge_energy[n_max_ge];
  Long64_t ge_ts[n_max_ge];

  Int_t n_spider = 0;
  Short_t spider_sid[n_max_spider];
  Short_t spider_ch[n_max_spider];
  Double_t spider_energy[n_max_spider];
  Long64_t spider_ts[n_max_spider];

  Int_t n_s3_sector = 0;
  Short_t s3_sector_sid[n_max_s3_sector];
  Short_t s3_sector_ch[n_max_s3_sector];
  Short_t s3_sector_id[n_max_s3_sector];
  Double_t s3_sector_energy[n_max_s3_sector];
  Long64_t s3_sector_ts[n_max_s3_sector];
  Int_t n_s3_ring = 0;
  Short_t s3_ring_sid[n_max_s3_ring];
  Short_t s3_ring_ch[n_max_s3_ring];
  Short_t s3_ring_id[n_max_s3_ring];
  Double_t s3_ring_energy[n_max_s3_ring];
  Long64_t s3_ring_ts[n_max_s3_ring];

  std::stringstream ss;
  ss.str("");
  ss << "random coincidence window [" << t_min_rand << ", " << t_max_rand << "] ns" ;
  TTree *tr = new TTree("tr_bg", ss.str().c_str());

  //
  memset(ge_sid, 0, sizeof(ge_sid));
  memset(ge_ch, 0, sizeof(ge_ch));
  memset(ge_energy, 0, sizeof(ge_energy));
  memset(ge_ts, 0, sizeof(ge_ts));

  memset(spider_sid, 0, sizeof(spider_sid));
  memset(spider_ch, 0, sizeof(spider_ch));
  memset(spider_energy, 0, sizeof(spider_energy));
  memset(spider_ts, 0, sizeof(spider_ts));

  memset(s3_sector_sid, 0, sizeof(s3_sector_sid));
  memset(s3_sector_ch, 0, sizeof(s3_sector_ch));
  memset(s3_sector_id, 0, sizeof(s3_sector_id));
  memset(s3_sector_energy, 0, sizeof(s3_sector_energy));
  memset(s3_sector_ts, 0, sizeof(s3_sector_ts));

  memset(s3_ring_sid, 0, sizeof(s3_ring_sid));
  memset(s3_ring_ch, 0, sizeof(s3_ring_ch));
  memset(s3_ring_id, 0, sizeof(s3_ring_id));
  memset(s3_ring_energy, 0, sizeof(s3_ring_energy));
  memset(s3_ring_ts, 0, sizeof(s3_ring_ts));

  //
  tr->Branch("n_ge", &n_ge, "n_ge/I");
  tr->Branch("ge_sid", ge_sid, "ge_sid[n_ge]/S");
  tr->Branch("ge_ch", ge_ch, "ge_ch[n_ge]/S");
  tr->Branch("ge_energy", ge_energy, "ge_energy[n_ge]/D");
  tr->Branch("ge_ts", ge_ts, "ge_ts[n_ge]/L");

  tr->Branch("n_spider", &n_spider, "n_spider/I");
  tr->Branch("spider_sid", spider_sid, "spider_sid[n_spider]/S");
  tr->Branch("spider_ch", spider_ch, "spider_ch[n_spider]/S");
  tr->Branch("spider_energy", spider_energy, "spider_energy[n_spider]/D");
  tr->Branch("spider_ts", spider_ts, "spider_ts[n_spider]/L");
  
  tr->Branch("n_s3_sector", &n_s3_sector, "n_s3_sector/I");
  tr->Branch("s3_sector_sid", s3_sector_sid, "s3_sector_sid[n_s3_sector]/S");
  tr->Branch("s3_sector_ch", s3_sector_ch, "s3_sector_ch[n_s3_sector]/S");
  tr->Branch("s3_sector_id", s3_sector_id, "s3_sector_id[n_s3_sector]/S");
  tr->Branch("s3_sector_energy", s3_sector_energy, "s3_sector_energy[n_s3_sector]/D");
  tr->Branch("s3_sector_ts", s3_sector_ts, "s3_sector_ts[n_s3_sector]/L");

  tr->Branch("n_s3_ring", &n_s3_ring, "n_s3_ring/I");
  tr->Branch("s3_ring_sid", s3_ring_sid, "s3_ring_sid[n_s3_ring]/S");
  tr->Branch("s3_ring_ch", s3_ring_ch, "s3_ring_ch[n_s3_ring]/S");
  tr->Branch("s3_ring_id", s3_ring_id, "s3_ring_id[n_s3_ring]/S");
  tr->Branch("s3_ring_energy", s3_ring_energy, "s3_ring_energy[n_s3_ring]/D");
  tr->Branch("s3_ring_ts", s3_ring_ts, "s3_ring_ts[n_s3_ring]/L");
  
  //
  Long64_t i_start = 0;

  Short_t sid1 = 0;
  Short_t ch1 = 0;
  Double_t energy1 = 0;
  Long64_t ts1 = 0;
  
  Short_t sid2 = 0;
  Short_t ch2 = 0;
  Double_t energy2 = 0;
  Long64_t ts2 = 0;

  while(true){//get first Ge data
    tr_in->GetEntry(i_start);
    energy = GetEnergy();
    ts_ns = GetTsns();
    i_start++;
    if(cid==0 && energy>cut_ge){
      break;
    }
  }
  
  Long64_t n_evt = 0;
  Long64_t i = i_start;
  Long64_t i_current = i-1;
  while(true){
    if(i==tr_in->GetEntries()) break;

    tr_in->GetEntry(i_current);
    energy = GetEnergy();
    ts_ns = GetTsns();

    sid1 = sid;
    ch1 = ch;
    energy1 = energy;
    ts1 = ts_ns;

    ge_sid[n_ge] = sid1;
    ge_ch[n_ge] = ch1; 
    ge_energy[n_ge] = energy1; 
    ge_ts[n_ge] = ts1; 
    n_ge++;

    // std::cout << std::endl;
    // std::cout << "info ge " << cid << " " << sid << " " << ch << " " << energy << " " << ts << std::endl;

    while(true){//search backward
      if(i_current--==0) break;
      tr_in->GetEntry(i_current);
      energy = GetEnergy();
      ts_ns = GetTsns();
      // std::cout << "info backward " << cid << " " << sid << " " << ch << " " << energy << " " << ts << std::endl;
      // std::cout << "time diff " << abs(ts-ts1) << std::endl;

      if(cid==0 && energy<=cut_ge){//if Ge data and small energy
        continue;
      }
      if(cid==1 && energy<=cut_si){//if Si/S3 data and small energy
        continue;
      }

      sid2 = sid;
      ch2 = ch;
      energy2 = energy;
      ts2 = ts_ns;
    
      if(abs(ts2-ts1)<=t_prompt){//for ge coincidence
        if(cid==0){
          ge_sid[n_ge] = sid2;
          ge_ch[n_ge] = ch2;
          ge_energy[n_ge] = energy2;
          ge_ts[n_ge] = ts2;
          n_ge++;
        }else{
          continue;
        }
      }
      else if(abs(ts2-ts1)>t_prompt && abs(ts2-ts1)<t_min_rand) continue;
      else if((abs(ts2-ts1))>=t_min_rand && abs(ts2-ts1)<=t_max_rand){//this coincidence
        if(cid==1){
          if(!(std::find(v_s3_sid.begin(),v_s3_sid.end(),sid) != v_s3_sid.end())){// spider data
            spider_sid[n_spider] = sid2;
            spider_ch[n_spider] = ch2;
            spider_energy[n_spider] = energy2;
            spider_ts[n_spider] = ts2;
            n_spider++;
          }else{
            if(sid>=11&&sid<=12){
              s3_ring_sid[n_s3_ring] = sid2;
              s3_ring_ch[n_s3_ring] = ch2;
              s3_ring_id[n_s3_ring] = map_s3_ring_id[10000*cid+100*sid+ch];
              s3_ring_energy[n_s3_ring] = energy2;
              s3_ring_ts[n_s3_ring] = ts2;
              n_s3_ring++;
            }else{
              s3_sector_sid[n_s3_sector] = sid2;
              s3_sector_ch[n_s3_sector] = ch2;
              s3_sector_id[n_s3_sector] = map_s3_sector_id[10000*cid+100*sid+ch];
              s3_sector_energy[n_s3_sector] = energy2;
              s3_sector_ts[n_s3_sector] = ts2;
              n_s3_sector++;
            }
          }
        }else continue;
      }else break;
    }//while search backward

    while(true){//search forward
      if(i>=tr_in->GetEntries()) break;
      tr_in->GetEntry(i);
      energy = GetEnergy();
      ts_ns = GetTsns();
      // std::cout << "info forward " << cid << " " << sid << " " << ch << " " << energy << " " << ts << std::endl;
      // std::cout << "time diff " << abs(ts-ts1) << std::endl;

      if(cid==0 && energy<=cut_ge){//if ge data and small energy
        i++;
        continue;
      }
      if(cid==1 && energy<=cut_si){//if si data and small energy
        i++;
        continue;
      }

      sid2 = sid;
      ch2 = ch;
      energy2 = energy;
      ts2 = ts_ns;

      if((abs(ts2-ts1))<=t_prompt){//for ge coincidence
        i++;
        if(cid==0){
          ge_sid[n_ge] = sid2;
          ge_ch[n_ge] = ch2;
          ge_energy[n_ge] = energy2;
          ge_ts[n_ge] = ts2;
          n_ge++;
        }else{
          continue;
        }
      }
      else if(abs(ts2-ts1)>t_prompt && abs(ts2-ts1)<t_min_rand){
        i++;
        continue;
      }
      else if((abs(ts2-ts1))>=t_min_rand && abs(ts2-ts1)<=t_max_rand){//this coincidence      
        i++;
        if(cid==1){
          if(!(std::find(v_s3_sid.begin(),v_s3_sid.end(),sid) != v_s3_sid.end())){// spider data
            spider_sid[n_spider] = sid2;
            spider_ch[n_spider] = ch2;
            spider_energy[n_spider] = energy2;
            spider_ts[n_spider] = ts2;
            n_spider++;
          }else{
            if(sid>=11&&sid<=12){
              s3_ring_sid[n_s3_ring] = sid2;
              s3_ring_ch[n_s3_ring] = ch2;
              s3_ring_id[n_s3_ring] = map_s3_ring_id[10000*cid+100*sid+ch];
              s3_ring_energy[n_s3_ring] = energy2;
              s3_ring_ts[n_s3_ring] = ts2;
              n_s3_ring++;
            }else{
              s3_sector_sid[n_s3_sector] = sid2;
              s3_sector_ch[n_s3_sector] = ch2;
              s3_sector_id[n_s3_sector] = map_s3_sector_id[10000*cid+100*sid+ch];
              s3_sector_energy[n_s3_sector] = energy2;
              s3_sector_ts[n_s3_sector] = ts2;
              n_s3_sector++;
            }
          }
        }else continue;
      }else break;
    } // while search forward

    // std::cout << "n_ge " << n_ge << std::endl;
    // std::cout << " n_spider " << n_spider << std::endl;
    // std::cout << " n_s3_sector " << n_s3_sector << " n_s3_ring " << n_s3_ring << std::endl;
    if(n_ge*(n_spider+n_s3_sector+n_s3_ring) > 0){
      n_evt++;

      file_out->cd();
      tr->Fill();

      if(n_evt%10000==0){
        std::cout << "\r" << n_evt << "  " << i << "/" << tr_in->GetEntries();
        std::cout << std::flush;
      }
    }

    n_ge = 0;
    memset(ge_sid, 0, sizeof(ge_sid));
    memset(ge_ch, 0, sizeof(ge_ch));
    memset(ge_energy, 0, sizeof(ge_energy));
    memset(ge_ts, 0, sizeof(ge_ts));

    n_spider = 0;
    memset(spider_sid, 0, sizeof(spider_sid));
    memset(spider_ch, 0, sizeof(spider_ch));
    memset(spider_energy, 0, sizeof(spider_energy));
    memset(spider_ts, 0, sizeof(spider_ts));

    n_s3_sector = 0;
    memset(s3_sector_sid, 0, sizeof(s3_sector_sid));
    memset(s3_sector_ch, 0, sizeof(s3_sector_ch));
    memset(s3_sector_id, 0, sizeof(s3_sector_id));
    memset(s3_sector_energy, 0, sizeof(s3_sector_energy));
    memset(s3_sector_ts, 0, sizeof(s3_sector_ts));

    n_s3_ring = 0;
    memset(s3_ring_sid, 0, sizeof(s3_ring_sid));
    memset(s3_ring_ch, 0, sizeof(s3_ring_ch));
    memset(s3_ring_id, 0, sizeof(s3_ring_id));
    memset(s3_ring_energy, 0, sizeof(s3_ring_energy));
    memset(s3_ring_ts, 0, sizeof(s3_ring_ts));

    while(true){//get next Ge data
      if(i>=tr_in->GetEntries()) break;
      tr_in->GetEntry(i);
      energy = GetEnergy();
      ts_ns = GetTsns();
      i++;
      if(cid==0 && energy>cut_ge){
        break;
      }
    }
    i_current = i-1;

  }//while

  std::cout << std::endl;

  file_out->cd();
  tr->Write();
}
//
double build::GetEnergy()
{
  int key = 10000*cid+100*sid+ch;

  auto it = map_cali_data.find(key);
  if(it==map_cali_data.end()){
    return 0;
  }

  double energy = 0;
  if(cid==0){
    timing_ge_energy_keV = (map_cali_data[key][0]+evte*map_cali_data[key][1]+evte*evte*map_cali_data[key][2])*map_k_data[key];
    energy = map_cali_data[key][0]+evte*map_cali_data[key][1]+evte*evte*map_cali_data[key][2]+rndm->Uniform(-0.5,0.5);
    energy *= map_k_data[key];
  }

  if(cid==1){
    energy = map_cali_data[key][0]+evte*map_cali_data[key][1]+evte*evte*map_cali_data[key][2]+rndm->Uniform(-5.,5.);

    if(map_s3_sector_id.find(key)!=map_s3_sector_id.end()){
      energy = map_s3_sector_cor_data[map_s3_sector_id[key]][0]+energy*map_s3_sector_cor_data[map_s3_sector_id[key]][1];
    }
    if(map_s3_ring_id.find(key)!=map_s3_ring_id.end()){
      energy = map_s3_ring_cor_data[map_s3_ring_id[key]][0]+energy*map_s3_ring_cor_data[map_s3_ring_id[key]][1];
    }
  }

  return energy;
}

//
Long64_t build::GetTsns()
{
  int key = 10000*cid+100*sid+ch;

  Long64_t ts_ns = 0;
  if(cid==0){
    if(cfdft){
      ts_ns = ts*8;
    }else{
      ts_ns = (ts*2-cfds+cfd/16384.)*4.;
    }

    if(time_alignment.total()){
      ts_ns = static_cast<Long64_t>(static_cast<double>(ts_ns)+time_alignment.TotalGeOffset(key, cfdft, timing_ge_energy_keV));
    }else{
    ts_ns += map_ts_offset[key];
    if(time_alignment.active()){
      double delta = 0.;
      if(time_alignment.GeDelta(key, cfdft, timing_ge_energy_keV, delta))
        ts_ns = static_cast<Long64_t>(static_cast<double>(ts_ns)+delta);
    }else{
      for(int i=0;i<tf->GetNpar();++i) tf->SetParameter(i, map_ge_tw_data[key][i]);
      ts_ns += tf->Eval(energy);
    }
    }
  }

  if(cid==1){
    if(sr==250) ts_ns = ts*8;
    if(sr==100) ts_ns = ts*10; 

    if(time_alignment.total()){
      ts_ns = static_cast<Long64_t>(static_cast<double>(ts_ns)+time_alignment.TotalParticleOffsetFromAdc(key, evte));
    }else{
      ts_ns += map_ts_offset[key];
      const bool is_s3 = map_s3_ring_id.count(key) || map_s3_sector_id.count(key);
      ts_ns = static_cast<Long64_t>(static_cast<double>(ts_ns)+time_alignment.ParticleDeltaFromAdc(key, evte, is_s3));
    }
    ts_ns += rndm->Uniform(-5,5);
  }

  return ts_ns;
}

//
void build::SaveFile()
{
  file_out->cd();
  file_out->Close();
}

//
bool build::ReadCaliData()
{
  std::cout << "start read cali data" << std::endl;

  //
  std::ifstream fi_cali_ge;
  fi_cali_ge.open("../pars/cali_ge.txt");
  if(!fi_cali_ge){
    std::cout << "can not open cali.txt" << std::endl;
    return 0;
  }

  std::string line;
  std::getline(fi_cali_ge, line);

  int cid, sid, ch;
  double par0, par1, par2, chi2;
  int key = 0;

  while(1){
    fi_cali_ge >> cid >> sid >> ch >> par0 >> par1 >> par2 >> chi2;
    if(!fi_cali_ge.good()) break;

    key = 10000*cid+100*sid+ch; // ch from 0 to 15
    std::vector<double> value;
    value.push_back(par0);
    value.push_back(par1);
    value.push_back(par2);
    value.push_back(chi2);

    map_cali_data.insert(std::pair<int, std::vector<double>>(key, value));
  }

  fi_cali_ge.close();

  //
  std::ifstream fi_k;
  fi_k.open(TString::Format("../pars/run_k/k_%04d.txt", run).Data());
  if(!fi_k){
    std::cout << "can not open k data" << std::endl;
    for(int i=2;i<=5;i++){
      for(int j=0;j<16;j++){
        map_k_data.insert(std::pair<int, double>(100*i+j, 1.));
      }
    }
  }else{
    std::string line;
    std::getline(fi_k, line);

    double k;

    while(1){
      fi_k >> cid >> sid >> ch >> k;
      if(!fi_k.good()) break;

      key = 10000*cid+100*sid+ch; // ch from 0 to 15

      map_k_data.insert(std::pair<int, double>(key, k));
    }

    fi_k.close();
  }

  //
  std::ifstream fi_cali_si;
  if(run==460 || (run>=472 && run<=605)){
    fi_cali_si.open("../pars/cali_si_r460.txt");
  }else if((run>=616 && run<=674) || run==722){
    fi_cali_si.open("../pars/cali_si_r457_low_current.txt");
  }else if(run==457 || (run>=680 && run<=718)){
    fi_cali_si.open("../pars/cali_si_r457_high_current.txt");
  }else{
    std::cout << "wrong run number." << std::endl;
    return 0;
  }

  if(!fi_cali_si){
    std::cout << "can not open si cali.txt" << std::endl;
    return 0;
  }

  std::getline(fi_cali_si, line);

  while(1){
    fi_cali_si >> cid >> sid >> ch >> par0 >> par1 >> par2 >> chi2;
    if(!fi_cali_si.good()) break;

    key = 10000*cid+100*sid+ch; // ch from 0 to 15
    std::vector<double> value;
    value.push_back(par0);
    value.push_back(par1);
    value.push_back(par2);
    value.push_back(chi2);

    map_cali_data.insert(std::pair<int, std::vector<double>>(key, value));
  }

  fi_cali_si.close();


  return 1;
}

//
void build::PrintCaliData()
{
  std::cout << "start print cali data" << std::endl;

  std::map<int, std::vector<double>>::iterator it_cali = map_cali_data.begin();
  for(;it_cali!=map_cali_data.end();it_cali++){
    std::cout << it_cali->first << " => " << it_cali->second[0] << " " << it_cali->second[1] << " " << it_cali->second[2] << '\n';
  }

  std::map<int, double>::iterator it_k = map_k_data.begin();
  for(;it_k!=map_k_data.end();it_k++){
    std::cout << it_k->first << " => " << it_k->second<< '\n';
  }
}

//
void build::ReadS3CorData() 
{
  int id;
  double par0, par1;

  std::string line;

  std::ifstream fi_fb_cor_sector(TString::Format("../pars/run_fb_cor/correction_sector_%04d.txt",run).Data());
  if(!fi_fb_cor_sector){
    std::cout << "can not open sector correction file." << std::endl;
    for(int i=1;i<=32;i++){
      map_s3_sector_cor_data[i] = {0, 1};
    }
  }else{
    std::getline(fi_fb_cor_sector, line);
    
    while(1){
      fi_fb_cor_sector >> id >> par0 >> par1;
      if(!fi_fb_cor_sector.good()) break;

      map_s3_sector_cor_data[id] = {par0, par1};
    }
    fi_fb_cor_sector.close();
  }

  std::ifstream fi_fb_cor_ring(TString::Format("../pars/run_fb_cor/correction_ring_%04d.txt",run).Data());
  if(!fi_fb_cor_ring){
    std::cout << "can not open ring correction file." << std::endl;
    for(int i=1;i<=24;i++){
      map_s3_ring_cor_data[i] = {0, 1};
    }
  }else{
    std::getline(fi_fb_cor_ring, line);
    
    while(1){
      fi_fb_cor_ring >> id >> par0 >> par1;
      if(!fi_fb_cor_ring.good()) break;

      map_s3_ring_cor_data[id] = {par0, par1};
    }
    fi_fb_cor_ring.close();
  }
}

//
void build::PrintS3CorData()
{
  std::cout << "\nstart print s3 cor data" << std::endl;

  std::cout << "\nfor sector ..." << std::endl;
  std::map<int, std::vector<double>>::iterator it_sector = map_s3_sector_cor_data.begin();
  for(;it_sector!=map_s3_sector_cor_data.end();it_sector++){
    std::cout << it_sector->first << " => " << it_sector->second[0] << " " << it_sector->second[1] << '\n';
  }

  std::cout << "\nfor ring ..." << std::endl;
  std::map<int, std::vector<double>>::iterator it_ring = map_s3_ring_cor_data.begin();
  for(;it_ring!=map_s3_ring_cor_data.end();it_ring++){
    std::cout << it_ring->first << " => " << it_ring->second[0] << " " << it_ring->second[1] << '\n';
  }
}

//
bool build::ReadTsOffset()
{
  std::cout << "start read ts offset data" << std::endl;

  std::ifstream fi;
  if((run>=457 && run<=462) || (run>=472 && run<=605)){
    fi.open("../pars/ts_offset_r472_r605.txt");
  }else if((run>=616 && run<=674) || (run>=721 && run<=723)){
    fi.open("../pars/ts_offset_r616_r674.txt");
  }else if(run>=680 && run<=718){
    fi.open("../pars/ts_offset_r680_r718.txt");
  }else{
    std::cout << "wrong run number." << std::endl;
    return 0;
  }
  
  if(!fi){
    std::cout << "can not open ts offset data" << std::endl;
    return 0;
  }

  std::string line;
  std::getline(fi, line);

  int cid, sid, ch;
  double ts_offset, chi2;
  int key = 0;

  while(1){
    fi >> cid >> sid >> ch >> ts_offset >> chi2;
    if(!fi.good()) break;

    key = 10000*cid+100*sid+ch; // ch from 0 to 15

    map_ts_offset.insert(std::pair<int, double>(key, ts_offset));
  }

  fi.close();

  return 1;
}

//
void build::PrintTsOffset()
{
  std::cout << "start print ts offset data" << std::endl;

  std::map<int, double>::iterator it = map_ts_offset.begin();
  for(;it!=map_ts_offset.end();it++){
    std::cout << it->first << " => " << it->second << '\n';
  }
}
  
//
bool build::ReadGeTwData()
{
  std::cout << "start read ge ts time walk data" << std::endl;

  std::ifstream fi;
  fi.open(TString::Format("../ge_time_walk/par/ge_time_walk_%s_%s.txt", TF, CUT).Data());

  if(!fi){
    std::cout << "can not open ge ts time walk data" << std::endl;
    return 0;
  }

  std::string line;
  std::getline(fi, line);

  int n_par = tf->GetNpar();
  std::cout << "n_par " << n_par << std::endl;

  int sid, ch;
  double par[n_par];
  int key = 0;

  while(1){
    fi >> sid >> ch;
    for(int i=0;i<n_par;++i) fi >> par[i];
    if(!fi.good()) break;

    key = 100*sid+ch; // ch from 0 to 15
    std::vector<double> value;
    for(int i=0;i<n_par;++i) value.push_back(par[i]);

    map_ge_tw_data.insert(std::pair<int, std::vector<double>>(key, value));
  }

  fi.close();

  return 1;
}

//
void build::PrintGeTwData()
{
  std::cout << "start print ge ts time walk data" << std::endl;

  std::map<int, std::vector<double>>::iterator it_ge_tw = map_ge_tw_data.begin();
  for(;it_ge_tw!=map_ge_tw_data.end();it_ge_tw++){
    std::cout << it_ge_tw->first << " => ";
    for(auto i=0u;i<it_ge_tw->second.size();++i) std::cout << it_ge_tw->second[i] << " ";
    std::cout << '\n';
  }
}

//
bool build::InitMapSectorRingID()
{
  std::cout << "start init ring and sector id." << std::endl;

  if((run>=457 && run<=462) || (run>=472 && run<=605) || (run>=680 && run<=718)){
    map_s3_sector_id = {
      {10900, 10},
      {10901, 12},
      {10902, 14},
      {10903, 16},
      {10904, 18},
      {10905, 20},
      {10906, 22},
      {10907, 24},
      {10908, 26},
      {10909, 28},
      {10910, 30},
      {10911, 32},
      {10912, 2},
      {10913, 4},
      {10914, 6},
      {10915, 8},
      {11000, 11},
      {11001, 13},
      {11002, 15},
      {11003, 17},
      {11004, 19},
      {11005, 21},
      {11006, 23},
      {11007, 25},
      {11008, 27},
      {11009, 29},
      {11010, 31},
      {11011, 1},
      {11012, 3},
      {11013, 5},
      {11014, 7},
      {11015, 9}
    };
  }else if((run>=616 && run<=674) || (run>=721 && run<=723)){
    map_s3_sector_id = {
      {10200, 10},
      {10201, 12},
      {10202, 14},
      {10203, 16},
      {10204, 18},
      {10205, 20},
      {10206, 22},
      {10207, 24},
      {10208, 26},
      {10209, 28},
      {10210, 30},
      {10211, 32},
      {10212, 2},
      {10213, 4},
      {10214, 6},
      {10215, 8},
      {10300, 11},
      {10301, 13},
      {10302, 15},
      {10303, 17},
      {10304, 19},
      {10305, 21},
      {10306, 23},
      {10307, 25},
      {10308, 27},
      {10309, 29},
      {10310, 31},
      {10311, 1},
      {10312, 3},
      {10313, 5},
      {10314, 7},
      {10315, 9}
    };
  }else{
  std::cout << "wrong run number." << std::endl;
    return 0;
  }

  map_s3_ring_id = {
    {11102, 2},
    {11103, 4},
    {11104, 6},
    {11105, 8},
    {11106, 10},
    {11107, 12},
    {11108, 14},
    {11109, 16},
    {11110, 18},
    {11111, 20},
    {11112, 22},
    {11113, 24},
    {11200, 1},
    {11201, 3},
    {11202, 5},
    {11203, 7},
    {11204, 9},
    {11205, 11},
    {11206, 13},
    {11207, 15},
    {11208, 17},
    {11209, 19},
    {11210, 21},
    {11211, 23}
  };

  return 1;
}

//
void build::PrintMapSectorRingID()
{
  std::cout << "start print ring and sector id." << std::endl;

  for(auto &it:map_s3_sector_id){
    std::cout << it.first << " => " << it.second << std::endl;
  }
  std::cout << std::endl;

  for(auto &it:map_s3_ring_id){
    std::cout << it.first << " => " << it.second << std::endl;
  }
}

