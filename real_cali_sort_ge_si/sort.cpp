#include "sort.h"

#include <vector>
#include "math.h"
#include "TNamed.h"

//
sort::sort(const std::string &filename_in, const std::string &filename_out, const int r)
{
  benchmark = new TBenchmark;

  run = r;

  file_in = TFile::Open(filename_in.c_str());
  if(file_in->IsZombie()){
    std::cout << "open file " << filename_in.c_str() << " error!" << std::endl;
    delete file_in;

    return ;
  }
  tr_in = (TTree*)file_in->Get("tree");
  rd = new raw(tr_in); 

  //
  if(std::strcmp(TF,"tf1") == 0) tf = new TF1("tf", "[0]+[1]/pow(x,0.5)+[2]/pow(x,1)+[3]/pow(x,2)", 0, 4096);
  else if(std::strcmp(TF, "tf2") == 0) tf = new TF1("tf","[0]+[1]*exp([2]/pow(x,1.5))", 0, 4096);
  else{
    throw std::invalid_argument("wrong tf name.");
  }

  //
  if(!InitMapSectorRingID()){
    throw std::invalid_argument("can not init ring and sector id data.");
  }
  PrintMapSectorRingID();

  // 
  if(!ReadCaliData()){
    throw std::invalid_argument("can not read cali data.");
  }
  PrintCaliData();
  
  if(!ReadS3AttData()){
    throw std::invalid_argument("can not read ts s3 att data.");
  }
  PrintS3AttData();

  if(!ReadS3CorData()){
    throw std::invalid_argument("can not read ts s3 cor data.");
  }
  PrintS3CorData();

  time_alignment.Load("../pars", run);
  // Complete profiles are self-contained: no external legacy offset is read.
  if(!time_alignment.total()){
    if(!ReadTSOffset()){
      throw std::invalid_argument("can not read ts offset data.");
    }
    PrintTSOffset();
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
  
  //
  file_out = TFile::Open(filename_out.c_str(), "recreate");
  file_out->cd();
  TNamed("sm150_time_alignment_profile", time_alignment.profile().c_str()).Write();
  TNamed("sm150_time_alignment_source_sha256", time_alignment.source_sha256().c_str()).Write();

  tr_out = new TTree("tree", "GDDAQ cali and sort data");
  tr_out->Branch("cid", &sd.cid);
  tr_out->Branch("sid", &sd.sid);
  tr_out->Branch("ch", &sd.ch);
  tr_out->Branch("evte", &sd.evte);
  tr_out->Branch("ts", &sd.ts);
}

//
sort::~sort()
{
  file_in->Close();
  file_out->Close();
}

//
void sort::Process()
{
  std::cout << "process " << std::endl;

  benchmark->Start("process");

  rd->GetEntry(rd->GetEntries()-1);
  sd(rd->cid, rd->sid, rd->ch, 0., 0);
  GetEnergy();
  if(map_cali.find(10000*rd->cid+100*rd->sid+rd->ch) == map_cali.end()){
    std::cout << "can not find cali data for this channel." << std::endl;
    std::cout << "cid " << rd->cid << " sid " << rd->sid << " ch " << rd->ch << std::endl;
    return;
  }

  GetTsns();
  std::cout << "total time " << sd.ts/1000000000 << "s" << std::endl; 

  Long64_t ts_save = TIMEBUFFER;
  Long64_t key;
  for(Long64_t i=0;i<rd->GetEntries();++i){
    rd->GetEntry(i);
    sd(rd->cid, rd->sid, rd->ch, 0., 0);

    if(map_cali.find(10000*rd->cid+100*rd->sid+rd->ch) == map_cali.end()){
      std::cout << "can not find cali data for this channel." << std::endl;
      std::cout << "cid " << rd->cid << " sid " << rd->sid << " ch " << rd->ch << std::endl;
      return;
    }
    
    GetEnergy();
    if(sd.evte<CUTGE) continue;
    GetTsns();

    key = ((((((sd.ts)<<4)+sd.cid)<<4)+sd.sid)<<6)+sd.ch;
    map_sort_data.insert(std::make_pair(key, sd));

    if(sd.ts>ts_save){
      Long64_t length = map_sort_data.size() * 0.9;
      for(Long64_t j=0;j<length;++j){
        it_sd = map_sort_data.begin();
        sd = it_sd->second;

        file_out->cd();
        tr_out->Fill();

        map_sort_data.erase(it_sd);
      }
    
      ts_save += TIMEBUFFER;
      std::cout << "\r" << sd.ts/1000000000 << "\r";
      std::cout << std::flush;
    }
  }

  int l = map_sort_data.size();
  for(int i=0;i<l;++i){
    it_sd = map_sort_data.begin();
    sd = it_sd->second;

    file_out->cd();
    tr_out->Fill();

    map_sort_data.erase(it_sd);
  }
  std::cout << "\r" << sd.ts/1000000000 << "\n";

  file_out->cd();
  tr_out->Write();
  file_out->Close();
  
  benchmark->Show("process");
}

//
void sort::GetTsns()
{
  int key = 10000*rd->cid+100*rd->sid+rd->ch;

  if(rd->cid==0){
    if(rd->cfdft){
      sd.ts = rd->ts*8;
    }else{
      sd.ts = (rd->ts*2-rd->cfds+rd->cfd/16384.)*4.;
    }
    
    if(time_alignment.total()){
      sd.ts = static_cast<Long64_t>(static_cast<double>(sd.ts)+time_alignment.TotalGeOffset(key, rd->cfdft, timing_ge_energy_keV));
    }else{
    sd.ts += map_ts_offset[key];
    if(time_alignment.active()){
      double delta = 0.;
      if(time_alignment.GeDelta(key, rd->cfdft, timing_ge_energy_keV, delta))
        sd.ts = static_cast<Long64_t>(static_cast<double>(sd.ts)+delta);
    }else{
#ifdef TW
      for(int i=0;i<tf->GetNpar();++i) tf->SetParameter(i, map_ge_tw[key][i]);
      sd.ts += tf->Eval(sd.evte);
#endif
    }
    }
  }

  if(rd->cid==1){
    if(rd->sr==250) sd.ts = rd->ts*8;
    if(rd->sr==100) sd.ts = rd->ts*10; 

    if(time_alignment.total()){
      sd.ts = static_cast<Long64_t>(static_cast<double>(sd.ts)+time_alignment.TotalParticleOffsetFromAdc(key, rd->evte));
    }else{
      sd.ts += map_ts_offset[key];
      const bool is_s3 = map_s3_ring_id.count(key) || map_s3_sector_id.count(key);
      sd.ts = static_cast<Long64_t>(static_cast<double>(sd.ts)+time_alignment.ParticleDeltaFromAdc(key, rd->evte, is_s3));
    }
    sd.ts += rndm->Uniform(-5,5);
  }
}

//
void sort::GetEnergy()
{
  int key = 10000*rd->cid+100*rd->sid+rd->ch;
  // Square ADC values in floating point to avoid signed-int overflow above 46340.
  const double adc = rd->evte;

  if(rd->cid==0){
    // The fitted lookup uses calibrated physical energy, without display dither.
    timing_ge_energy_keV = (map_cali[key][0]+adc*map_cali[key][1]+adc*adc*map_cali[key][2])*map_k[key];
    sd.evte = map_cali[key][0]+adc*map_cali[key][1]+adc*adc*map_cali[key][2]+rndm->Uniform(-0.5,0.5);
    sd.evte *= map_k[key];
  }

  if(rd->cid==1){
    sd.evte = map_cali[key][0]+adc*map_cali[key][1]+adc*adc*map_cali[key][2]+rndm->Uniform(-5.,5.);

    // modified
    if(rd->sid==11) sd.evte *= 10.;
    if((run>475&&run<616) || (run>=680&&run<=718)){
      if(rd->sid==2){
        sd.evte *= 4.;
      } 
    }
    if((run>=616&&run<=674) || (run>=721 && run<=723)){
      if(rd->sid==9){
        sd.evte *= 4.;
      } 
    }

    // s3att
    if(map_s3_sector_id.find(key)!=map_s3_sector_id.end()){
      sd.evte *= map_s3_sector_att[map_s3_sector_id[key]];
    }
    if(map_s3_ring_id.find(key)!=map_s3_ring_id.end()){
      sd.evte *= map_s3_ring_att[map_s3_ring_id[key]];
    }
 
    // s3cor
    if(map_s3_sector_id.find(key)!=map_s3_sector_id.end()){
      sd.evte = map_s3_sector_cor[map_s3_sector_id[key]][0]+sd.evte*map_s3_sector_cor[map_s3_sector_id[key]][1];
    }
    if(map_s3_ring_id.find(key)!=map_s3_ring_id.end()){
      sd.evte = map_s3_ring_cor[map_s3_ring_id[key]][0]+sd.evte*map_s3_ring_cor[map_s3_ring_id[key]][1];
    }
  }
}

//
bool sort::ReadCaliData()
{
  std::cout << "start read cali data" << std::endl;

  //
  std::ifstream fi_cali_ge;
  if(run==722) fi_cali_ge.open("../cali_ge/cali_0722.txt");
  else fi_cali_ge.open("../pars/cali_ge.txt");
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

    map_cali.insert(std::pair<int, std::vector<double>>(key, value));
  }

  fi_cali_ge.close();

  //
  std::ifstream fi_k;
  fi_k.open(TString::Format("../pars/run_k/k_%04d.txt", run).Data());
  if(!fi_k){
    std::cout << "can not open k data" << std::endl;
    for(int i=2;i<=5;i++){
      for(int j=0;j<16;j++){
        map_k.insert(std::pair<int, double>(100*i+j, 1.));
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

      map_k.insert(std::pair<int, double>(key, k));
    }

    fi_k.close();
  }

  //
  std::ifstream fi_cali_si;
  if(run==377 || run==460 || (run>=472 && run<=605)){
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

    map_cali.insert(std::pair<int, std::vector<double>>(key, value));
  }

  fi_cali_si.close();

  return 1;
}

//
void sort::PrintCaliData()
{
  std::cout << "start print cali data" << std::endl;

  std::map<int, std::vector<double>>::iterator it_cali = map_cali.begin();
  for(;it_cali!=map_cali.end();it_cali++){
    std::cout << it_cali->first << " => " << it_cali->second[0] << " " << it_cali->second[1] << " " << it_cali->second[2] << '\n';
  }

  std::map<int, double>::iterator it_k = map_k.begin();
  for(;it_k!=map_k.end();it_k++){
    std::cout << it_k->first << " => " << it_k->second<< '\n';
  }
}

//
bool sort::ReadS3CorData() 
{
  int id;
  double par0, par1;

  if(run>=721 && run<=723){
    for(int i=1;i<=32;i++) map_s3_sector_cor[i] = {0, 1.};
    for(int i=1;i<=24;i++) map_s3_ring_cor[i] = {0, 1.};

    return 1;
  }

  std::string line;

  std::ifstream fi_fb_cor_sector(TString::Format("../pars/run_fb_cor/correction_sector_%04d.txt",run).Data());
  if(!fi_fb_cor_sector){
    std::cout << "can not open sector correction file." << std::endl;
    return 0;
  }else{
    std::getline(fi_fb_cor_sector, line);
    
    while(1){
      fi_fb_cor_sector >> id >> par0 >> par1;
      if(!fi_fb_cor_sector.good()) break;

      map_s3_sector_cor[id] = {par0, par1};
    }
    fi_fb_cor_sector.close();
  }

  std::ifstream fi_fb_cor_ring(TString::Format("../pars/run_fb_cor/correction_ring_%04d.txt",run).Data());
  if(!fi_fb_cor_ring){
    std::cout << "can not open ring correction file." << std::endl;
    return 0;
  }else{
    std::getline(fi_fb_cor_ring, line);
    
    while(1){
      fi_fb_cor_ring >> id >> par0 >> par1;
      if(!fi_fb_cor_ring.good()) break;

      map_s3_ring_cor[id] = {par0, par1};
    }
    fi_fb_cor_ring.close();
  }

  return 1;
}

//
void sort::PrintS3CorData()
{
  std::cout << "\nstart print s3 cor data" << std::endl;

  std::cout << "\nfor sector ..." << std::endl;
  std::map<int, std::vector<double>>::iterator it_sector = map_s3_sector_cor.begin();
  for(;it_sector!=map_s3_sector_cor.end();it_sector++){
    std::cout << it_sector->first << " => " << it_sector->second[0] << " " << it_sector->second[1] << '\n';
  }

  std::cout << "\nfor ring ..." << std::endl;
  std::map<int, std::vector<double>>::iterator it_ring = map_s3_ring_cor.begin();
  for(;it_ring!=map_s3_ring_cor.end();it_ring++){
    std::cout << it_ring->first << " => " << it_ring->second[0] << " " << it_ring->second[1] << '\n';
  }
}

//
bool sort::ReadS3AttData()
{
  int r;
  double par;

  bool flag = 0;

  if(run>=721 && run<=723){
    for(int i=1;i<=32;i++) map_s3_sector_att[i] = 1.;
    for(int i=1;i<=24;i++) map_s3_ring_att[i] = 1.;

    return 1;
  }

  std::ifstream fi_sector;
  for(int i=0;i<32;i++){
    flag = 0;
    fi_sector.open(Form("../pars/run_s3_att/sector%02d.txt",i+1));
    if(!fi_sector){
      std::cout << "can not open sector att file." << std::endl;
      return 0;
    }else{
      while(fi_sector >> r >> par){
        if(r==run){
          map_s3_sector_att[i+1] = par;
          flag = 1;
          break;
        }else continue;
      }
      fi_sector.close();
    }
  }
  if(!flag){
    std::cout << "can not find run " << run << " s3 att data." << std::endl;
    return 0;
  }

  std::ifstream fi_ring;
  for(int i=0;i<24;i++){
    flag = 0;
    fi_ring.open(Form("../pars/run_s3_att/ring%02d.txt",i+1));
    if(!fi_ring){
      std::cout << "can not open ring att file." << std::endl;
      return 0;
    }else{
      while(fi_ring >> r >> par){
        if(r==run){
          map_s3_ring_att[i+1] = par;
          flag = 1;
          break;
        }else continue;
      }
      fi_ring.close();
    }
  }
  if(!flag){
    std::cout << "can not find run " << run << " s3 att data." << std::endl;
    return 0;
  }

  return 1;
}

//
void sort::PrintS3AttData()
{
  std::cout << "\nstart print s3 att data" << std::endl;

  std::cout << "\n run number " << run << " for sector ..." << std::endl;
  std::map<int, double>::iterator it_sector = map_s3_sector_att.begin();
  for(;it_sector!=map_s3_sector_att.end();it_sector++){
    std::cout << it_sector->first << " => " << it_sector->second << '\n';
  }

  std::cout << "\n run number " << run << " for ring ..." << std::endl;
  std::map<int, double>::iterator it_ring = map_s3_ring_att.begin();
  for(;it_ring!=map_s3_ring_att.end();it_ring++){
    std::cout << it_ring->first << " => " << it_ring->second << '\n';
  }
}

//
bool sort::ReadTSOffset()
{
  std::cout << "start read ts offset data" << std::endl;

  std::ifstream fi;
  if(run==377 || (run>=457 && run<=462) || (run>=472 && run<=605)){
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
void sort::PrintTSOffset()
{
  std::cout << "start print ts offset data" << std::endl;

  std::map<int, double>::iterator it = map_ts_offset.begin();
  for(;it!=map_ts_offset.end();it++){
    std::cout << it->first << " => " << it->second << '\n';
  }
}

//
bool sort::ReadGeTwData()
{
  std::cout << "start read ge ts time walk data" << std::endl;

  std::ifstream fi;
  fi.open(TString::Format("../pars/ge_time_walk_%s.txt", TF).Data());

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

    map_ge_tw.insert(std::pair<int, std::vector<double>>(key, value));
  }

  fi.close();

  return 1;
}

//
void sort::PrintGeTwData()
{
  std::cout << "start print ge ts time walk data" << std::endl;

  std::map<int, std::vector<double>>::iterator it_ge_tw = map_ge_tw.begin();
  for(;it_ge_tw!=map_ge_tw.end();it_ge_tw++){
    std::cout << it_ge_tw->first << " => ";
    for(auto i=0u;i<it_ge_tw->second.size();++i) std::cout << it_ge_tw->second[i] << " ";
    std::cout << '\n';
  }
}

//
bool sort::InitMapSectorRingID()
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
    {11102, 1},
    {11103, 3},
    {11104, 5},
    {11105, 7},
    {11106, 9},
    {11107, 11},
    {11108, 13},
    {11109, 15},
    {11110, 17},
    {11111, 19},
    {11112, 21},
    {11113, 23},
    {11200, 2},
    {11201, 4},
    {11202, 6},
    {11203, 8},
    {11204, 10},
    {11205, 12},
    {11206, 14},
    {11207, 16},
    {11208, 18},
    {11209, 20},
    {11210, 22},
    {11211, 24}
  };

  return 1;
}

//
void sort::PrintMapSectorRingID()
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

