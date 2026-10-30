// Channel singles spectra from decoded C1/C2 or explicitly calibrated trees.
// Reuse formal energy parameters and physical mappings; never apply timing or coincidence gates.
#include "sm150_event_calibration.h"
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TH1D.h>
#include <TObjString.h>
#include <TROOT.h>
#include <chrono>
#include <iostream>
#include <memory>
#include <limits>

using namespace sm150::diagnostic;
struct Options {
  int run=0, segments=20;
  Long64_t sample=0;
  fs::path project,input,output;
  std::string mode="both",format="decoded",detector="all";
};
static Long64_t integer(const std::string& s) {
  require(!s.empty() && s.find_first_not_of("0123456789")==std::string::npos,"Expected a nonnegative integer: "+s);
  size_t end=0; const auto n=std::stoll(s,&end);require(end==s.size(),"Invalid integer");return n;
}
static Options options(int argc,char** argv) {
  Options o;
  for(int i=1;i<argc;++i) {
    const std::string k=argv[i];
    if(k=="--help") {
      std::cout<<"Usage: histogram --run RUN --output NEW.root [--input C2.root] [--project DIR]\n"
        <<"  --sample-entries N (0: full input, default) --segments N (default: 20)\n"
        <<"  --mode both|raw|calibrated --detector all|ge|particle|si|s3\n"
        <<"  --format decoded|calibrated (explicit Double_t input; no recalibration)\n"
        <<"Decoded input: tree/cid,sid,ch Short_t; evte UShort_t; pileup,outofr Bool_t.\n"
        <<"Calibrated input: tree/cid,sid,ch Short_t; evte Double_t in keV; quality unavailable.\n"
        <<"Outputs refuse overwrite. Full physical mapping/calibration supports the three in-beam periods.\n";
      std::exit(0);
    }
    require(i+1<argc,"Missing value for "+k);const std::string v=argv[++i];
    if(k=="--run") { auto n=integer(v);require(n>0 && n<10000,"Invalid run");o.run=int(n); }
    else if(k=="--input")o.input=v;
    else if(k=="--output")o.output=v;
    else if(k=="--project")o.project=v;
    else if(k=="--sample-entries")o.sample=integer(v);
    else if(k=="--segments") {auto n=integer(v);require(n>0 && n<=10000,"Invalid segments");o.segments=int(n);}
    else if(k=="--mode")o.mode=v;
    else if(k=="--format")o.format=v;
    else if(k=="--detector")o.detector=v;
    else throw std::runtime_error("Unknown option: "+k);
  }
  require(o.run>0 && !o.output.empty(),"--run and --output are required");
  require(o.mode=="both" || o.mode=="raw" || o.mode=="calibrated","Invalid mode");
  require(o.format=="decoded" || o.format=="calibrated","Invalid format");
  require(o.detector=="all" || o.detector=="ge" || o.detector=="particle" || o.detector=="si" || o.detector=="s3","Invalid detector");
  require(o.format!="calibrated" || o.mode=="calibrated","Calibrated input requires --mode calibrated");
  return o;
}
static bool present(const fs::path& p) {return fs::exists(p)||fs::is_symlink(p);}
static std::string key_name(int key) {std::ostringstream s;s<<"ch"<<std::setw(5)<<std::setfill('0')<<key;return s.str();}
static const char* kind_name(int kind) {return kind==1?"Ge":kind==2?"Si array":kind==3?"S3 ring":kind==4?"S3 sector":"unmapped";}
struct Spectra {
  int kind=0;Long64_t seen=0,clean=0,rejected=0,pileup=0,outofr=0,saturated=0;
  std::map<std::string,std::unique_ptr<TH1D>> h;
  void add(const std::string& name,int n,double lo,double hi,const std::string& axis) {
    auto p=std::make_unique<TH1D>(name.c_str(),(";"+axis+";Counts / bin").c_str(),n,lo,hi);
    p->SetDirectory(nullptr);h.emplace(name,std::move(p));
  }
  Spectra(int k,bool raw,bool energy,bool decoded):kind(k) {
    if(raw) {add("adc_all",32768,0,65536,"ADC channel");add("adc_clean",32768,0,65536,"ADC channel");}
    if(energy && kind) {
      const int n=kind==1?4096:2000;const double high=kind==1?4096:100000;
      if(decoded) {
        add("energy_base",n,0,high,"Energy (keV)");
        if(kind>=3)add("energy_pre_fb",n,0,high,"Energy (keV)");
      }
      add("energy_final",n,0,high,"Energy (keV)");
    }
  }
};
int main(int argc,char** argv) {
  try {
    const auto begin=std::chrono::steady_clock::now();const Options o=options(argc,argv);
    const auto project=fs::canonical(o.project.empty()?fs::canonical("/proc/self/exe").parent_path().parent_path():o.project);
    const auto input=fs::canonical(o.input.empty()?fs::path("/mnt/e/Sm150/rootfile/data_C2_"+rn(o.run)+"_"+rn(o.run)+".root"):o.input);
    const auto output=fs::absolute(o.output).lexically_normal();const fs::path partial=output.string()+".partial",receipt=output.string()+".json",receipt_partial=receipt.string()+".partial";
    for(const auto& p:{output,partial,receipt,receipt_partial})require(!present(p),"Refusing existing output: "+p.string());
    const bool decoded=o.format=="decoded", raw=o.mode!="calibrated", energy=o.mode!="raw";
    // Pure raw mode can also inspect source runs outside the in-beam periods.
    const bool inbeam=(o.run>=472&&o.run<=605)||(o.run>=616&&o.run<=674)||(o.run>=680&&o.run<=718);
    std::unique_ptr<Parameters> pars;
    if(inbeam)pars=std::make_unique<Parameters>(project,o.run,true);
    require(pars || (!energy && (o.detector=="all"||o.detector=="ge"||o.detector=="particle")),"Energy/physical detector modes require an in-beam run");
    const Stat source_before(input);std::unique_ptr<TFile> fi(TFile::Open(input.c_str(),"READ"));
    require(fi && !fi->IsZombie() && !fi->TestBit(TFile::kRecovered),"Cannot open a complete input ROOT");
    auto* tree=dynamic_cast<TTree*>(fi->Get("tree"));require(tree,"Input must contain decoded/calibrated tree; event build trees are not singles");
    require(tree->GetEntries()>0,"Empty input tree");
    std::map<std::string,std::string> types={{"cid","Short_t"},{"sid","Short_t"},{"ch","Short_t"},{"evte",decoded?"UShort_t":"Double_t"}};
    if(decoded){types["pileup"]="Bool_t";types["outofr"]="Bool_t";}
    tree->SetBranchStatus("*",false);
    for(const auto& kv:types) {
      auto* leaf=tree->GetLeaf(kv.first.c_str());require(leaf && kv.second==leaf->GetTypeName() && leaf->GetLenStatic()==1 && !leaf->GetLeafCount(),"Wrong scalar branch: "+kv.first+" (expected "+kv.second+")");
      tree->SetBranchStatus(kv.first.c_str(),true);
    }
    Short_t cid=0,sid=0,ch=0;UShort_t adc=0;Double_t already=0;Bool_t pileup=false,outofr=false;
    require(tree->SetBranchAddress("cid",&cid)>=0 && tree->SetBranchAddress("sid",&sid)>=0 && tree->SetBranchAddress("ch",&ch)>=0,"Cannot bind channel branches");
    require((decoded?tree->SetBranchAddress("evte",&adc):tree->SetBranchAddress("evte",&already))>=0,"Cannot bind energy");
    if(decoded)require(tree->SetBranchAddress("pileup",&pileup)>=0 && tree->SetBranchAddress("outofr",&outofr)>=0,"Cannot bind quality");
    tree->SetCacheSize(8*1024*1024);for(const auto& kv:types)tree->AddBranchToCache(kv.first.c_str(),true);tree->StopCacheLearningPhase();
    const Long64_t total=tree->GetEntries(),wanted=o.sample?std::min(o.sample,total):total;
    const int blocks=wanted==total?1:int(std::min<Long64_t>(o.segments,wanted));
    json segments=json::array();std::map<int,std::unique_ptr<Spectra>> spectra;
    auto selected=[&](int key,int kind) {
      return o.detector=="all"||(o.detector=="ge"&&key/10000==0)||(o.detector=="particle"&&key/10000==1)||
        (o.detector=="si"&&kind==2)||(o.detector=="s3"&&kind>=3);
    };
    auto make=[&](int key,int kind)->Spectra& {
      if(!spectra.count(key))spectra[key]=std::make_unique<Spectra>(kind,raw,energy,decoded);
      return *spectra.at(key);
    };
    if(pars)for(const auto& kv:pars->channels)if(selected(kv.first,kv.second.kind))make(kv.first,kv.second.kind);
    Long64_t processed=0,filtered=0,invalid=0,unmapped=0;Digest digest;std::vector<char> packed;packed.reserve(1<<20);
    auto pack=[&](const auto& v){const char* p=reinterpret_cast<const char*>(&v);packed.insert(packed.end(),p,p+sizeof(v));};
    Long64_t remaining=wanted;
    for(int b=0;b<blocks;++b) {
      const Long64_t lo=total/blocks*b+(total%blocks)*b/blocks,hi=total/blocks*(b+1)+(total%blocks)*(b+1)/blocks;
      const Long64_t left=blocks-b,n=std::min(hi-lo,remaining/left+(remaining%left?1:0));
      remaining-=n;const Long64_t first=lo+(hi-lo-n)/2,last=first+n;
      require(first>=lo && last<=hi,"Invalid sampling segment");segments.push_back({{"begin",first},{"end",last}});tree->SetCacheEntryRange(first,last);
      for(Long64_t i=first;i<last;++i) {
        require(tree->GetEntry(i)>0,"Read failure at entry "+std::to_string(i));++processed;
        pack(i);pack(cid);pack(sid);pack(ch);
        if(decoded){pack(adc);pack(pileup);pack(outofr);}else pack(already);
        if(packed.size()>=(1<<20)){digest.add(packed.data(),packed.size());packed.clear();}
        if(cid<0||cid>1||sid<0||sid>99||ch<0||ch>=16){++invalid;continue;}
        const int key=cid*10000+sid*100+ch;const Channel* c=nullptr;
        if(pars){auto it=pars->channels.find(key);if(it!=pars->channels.end())c=&it->second;}
        const int kind=c?c->kind:0;
        if(!selected(key,kind)){++filtered;continue;}
        auto& s=make(key,kind);++s.seen;if(!kind)++unmapped;
        if(raw)s.h.at("adc_all")->Fill(adc);
        const bool saturated=decoded&&adc==65535;
        s.pileup+=pileup;s.outofr+=outofr;s.saturated+=saturated;
        if(pileup||outofr||saturated){++s.rejected;continue;}
        ++s.clean;if(raw)s.h.at("adc_clean")->Fill(adc);
        if(!energy||!c)continue;
        if(!decoded){require(std::isfinite(already),"Nonfinite calibrated input");s.h.at("energy_final")->Fill(already);continue;}
        const double x=adc,base=(c->cal[0]+c->cal[1]*x+c->cal[2]*x*x)*c->gain;
        const double pre=c->energy(x,false),final=c->energy(x,true);
        require(std::isfinite(base)&&std::isfinite(pre)&&std::isfinite(final),"Nonfinite energy");
        s.h.at("energy_base")->Fill(base);s.h.at("energy_final")->Fill(final);
        if(kind>=3)s.h.at("energy_pre_fb")->Fill(pre);
      }
      std::cout<<"run "<<rn(o.run)<<" segment "<<b+1<<"/"<<blocks<<" entries "<<processed<<"/"<<wanted<<std::endl;
    }
    if(!packed.empty())digest.add(packed.data(),packed.size());
    require(processed==wanted,"Incomplete input");if(pars)pars->verify();require(source_before.same(Stat(input)),"Source changed during analysis");
    json meta={{"schema","sm150_raw2hist_v1"},{"status","complete"},{"run",rn(o.run)},
      {"period",pars?pars->group:"source"},{"input",input.string()},{"input_size_bytes",source_before.value.st_size},
      {"input_mtime_sec",source_before.value.st_mtim.tv_sec},{"input_mtime_nsec",source_before.value.st_mtim.tv_nsec},
      {"input_uuid",fi->GetUUID().AsString()},{"input_entries",total},{"entries_read",processed},{"sample_fraction",double(processed)/total},
      {"segments",segments},{"branches",types},{"sampled_values_sha256",digest.finish()},
      {"digest_layout",decoded?"native little-endian entry:int64,cid/sid/ch:int16,ADC:uint16,pileup/outofr:bool in entry order":"native little-endian entry:int64,cid/sid/ch:int16,energy:float64 in entry order"},
      {"mode",o.mode},{"format",o.format},{"detector",o.detector},{"invalid_address_entries",invalid},{"filtered_entries",filtered},{"unmapped_entries",unmapped},
      {"energy_mode",decoded?"formal_current_tables":"input_energy_no_recalibration"},
      {"selection",decoded?"adc_all includes quality flags; clean excludes pileup OR outofr OR ADC==65535; no time/energy/coincidence gate; normal and failed CFD combined":"quality unavailable in explicitly pre-calibrated input; no further selection"},
      {"energy_definitions",{{"base","(c0+c1*ADC+c2*ADC^2)*hardware_gain"},{"final","fb0+fb1*base*k_ge*attenuation; keV; no random jitter"}}},
      {"calibration_usage",decoded?"applied to decoded ADC":"tables used for mapping only; original input calibration version is unspecified"},
      {"parameter_snapshots",pars?pars->snapshots:json::object()},{"parameter_records",pars?pars->records:json::object()},
      {"root_version",gROOT->GetVersion()},{"binary_sha256",hash(read_text(fs::canonical("/proc/self/exe")))},
      {"source_code_sha256",hash(read_text(project/"raw2hist/main.cpp"))},
      {"shared_header_sha256",hash(read_text(project/"include/sm150_event_calibration.h"))},
      {"channels",json::array()}};
    for(const auto& kv:spectra) {
      const auto& s=*kv.second;json h=json::object();
      for(const auto& item:s.h){const auto& v=*item.second;h[item.first]={{"entries",v.GetEntries()},{"integral",v.Integral()},{"underflow",v.GetBinContent(0)},{"overflow",v.GetBinContent(v.GetNbinsX()+1)},{"bins",v.GetNbinsX()},{"low",v.GetXaxis()->GetXmin()},{"high",v.GetXaxis()->GetXmax()}};}
      meta["channels"].push_back({{"key",kv.first},{"directory",key_name(kv.first)},{"kind",s.kind},{"detector",kind_name(s.kind)},
        {"seen",s.seen},{"clean",s.clean},{"rejected",s.rejected},{"pileup",s.pileup},{"outofr",s.outofr},{"saturated",s.saturated},{"histograms",h}});
    }
    fs::create_directories(output.parent_path());std::unique_ptr<TFile> fo(TFile::Open(partial.c_str(),"CREATE","Channel singles spectra",404));
    require(fo && !fo->IsZombie(),"Cannot create output");
    for(const auto& kv:spectra){auto* dir=fo->mkdir(key_name(kv.first).c_str());require(dir,"Cannot make ROOT directory");dir->cd();for(const auto& item:kv.second->h)require(item.second->Write()>0,"Histogram write failed");}
    meta["elapsed_seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    fo->cd();require(TObjString(meta.dump().c_str()).Write("metadata")>0,"Metadata write failed");
    fo->Close();require(!fo->TestBit(TFile::kWriteError),"ROOT write error");
    {std::ofstream out(receipt_partial);out<<meta.dump(2)<<'\n';out.close();require(bool(out),"Receipt write failed");}
    // Reserve final names atomically without replacing an existing file.
    fs::create_hard_link(partial,output);fs::create_hard_link(receipt_partial,receipt);fs::remove(partial);fs::remove(receipt_partial);
    std::cout<<"Complete: "<<output<<"; "<<processed<<" entries, "<<meta["elapsed_seconds"]<<" s\n";
    return 0;
  } catch(const std::exception& e) {std::cerr<<"raw2hist: "<<e.what()<<'\n';return 1;}
}
