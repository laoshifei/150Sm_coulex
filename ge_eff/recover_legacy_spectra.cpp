// Reconstruct source spectra using explicitly documented historical assumptions.
// This is an isolated efficiency study, not a replacement for production build.
#include "sm150_event_calibration.h"
#include <TFile.h>
#include <TH1D.h>
#include <TLeaf.h>
#include <TObjString.h>
#include <TRandom3.h>
#include <TTree.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>
#include <tuple>

using namespace sm150::diagnostic;
struct LegacyHit {
  Long64_t time, entry;
  double energy;
  int key;
  bool clean;
};

static int selected_ring(int key) {
  static const std::map<int,int> groups = {
    {204,3},{205,3},{206,3},{212,3},{213,3},{214,3},{215,3},
    {211,4},{304,4},{305,4},{306,4},{308,4},{507,4},{508,4},{510,4},
    {209,5},{310,5},{311,5},{506,5},{509,5}};
  auto it=groups.find(key); return it==groups.end()?0:it->second;
}

int main(int argc, char** argv) {
  try {
    require(argc==4,"Usage: recover_legacy_spectra INPUT_C1.root CONFIG.json NEW_OUTPUT.root");
    fs::path input=fs::canonical(argv[1]), conf=fs::canonical(argv[2]);
    fs::path output=fs::absolute(argv[3]);
    require(!fs::exists(output)&&!fs::exists(output.string()+".json"),"Output already exists");
    const Stat original(input);
    const auto start=std::chrono::steady_clock::now();
    std::string text=read_text(conf);json config=json::parse(text);
    require(config.at("run")=="0722","Expected source run 0722");
    std::map<int,std::array<double,3>> coefficients;
    std::map<int,double> offsets;
    for(auto it=config.at("channels").begin();it!=config.at("channels").end();++it) {
      int key=std::stoi(it.key());
      coefficients[key]=it.value().at("coefficients").get<std::array<double,3>>();
      offsets[key]=it.value().at("time_offset_ns").get<double>();
    }
    std::unique_ptr<TFile> file(TFile::Open(input.c_str(),"READ"));
    require(file&&!file->IsZombie()&&!file->TestBit(TFile::kRecovered),"Invalid source file");
    auto* tree=dynamic_cast<TTree*>(file->Get("tree"));require(tree,"Missing tree");
    tree->SetBranchStatus("*",false);
    const std::map<std::string,std::string> types={{"cid","Short_t"},{"sid","Short_t"},{"ch","Short_t"},
      {"evte","UShort_t"},{"ts","Long64_t"},{"sr","Short_t"},{"cfd","Short_t"},{"cfds","Short_t"},
      {"cfdft","Bool_t"},{"pileup","Bool_t"},{"outofr","Bool_t"}};
    for(const auto& kv:types) {
      auto* leaf=tree->GetLeaf(kv.first.c_str());
      require(leaf&&leaf->GetTypeName()==kv.second&&!leaf->GetLeafCount()&&leaf->GetLenStatic()==1,
              "Unexpected branch type: "+kv.first);
      tree->SetBranchStatus(kv.first.c_str(),true);
    }
    Short_t cid=0,sid=0,ch=0,sr=0,cfd=0,cfds=0;
    UShort_t adc=0;Long64_t ts=0;Bool_t failed=false,pileup=false,outofr=false;
    tree->SetBranchAddress("cid",&cid);tree->SetBranchAddress("sid",&sid);
    tree->SetBranchAddress("ch",&ch);tree->SetBranchAddress("sr",&sr);
    tree->SetBranchAddress("evte",&adc);tree->SetBranchAddress("ts",&ts);
    tree->SetBranchAddress("cfd",&cfd);tree->SetBranchAddress("cfds",&cfds);
    tree->SetBranchAddress("cfdft",&failed);tree->SetBranchAddress("pileup",&pileup);
    tree->SetBranchAddress("outofr",&outofr);
    tree->SetCacheSize(32*1024*1024);
    for(const auto& kv:types)tree->AddBranchToCache(kv.first.c_str(),true);
    tree->StopCacheLearningPhase();
    std::map<std::string,std::unique_ptr<TH1D>> hist;
    auto make=[&](const std::string& name) {
      auto h=std::make_unique<TH1D>(name.c_str(),";Energy (keV);Counts / 1 keV",4096,0,4096);
      h->SetDirectory(nullptr);hist[name]=std::move(h);
    };
    for(const auto& kv:coefficients)if(selected_ring(kv.first)) {
      make("ge"+std::to_string(kv.first)+"_all");make("ge"+std::to_string(kv.first)+"_clean");
    }
    std::vector<LegacyHit> hits;hits.reserve(size_t(tree->GetEntries()));
    TRandom3 random(config.at("dither_seed").get<unsigned>());
    json counts=json::object();Long64_t previous=-1;
    for(Long64_t i=0;i<tree->GetEntries();++i) {
      require(tree->GetEntry(i)>0,"Input read failure");
      require(cid==0&&sr==250,"Expected Ge-only 250-MHz source data");
      require(ts>=previous,"Coarse timestamps decreased");previous=ts;
      int key=100*sid+ch;require(coefficients.count(key),"Missing calibration for "+std::to_string(key));
      const auto& p=coefficients.at(key);
      double x=adc,e=p[0]+p[1]*x+p[2]*x*x+random.Uniform(-.5,.5);
      require(std::isfinite(e),"Nonfinite energy");
      bool clean=!pileup&&!outofr&&adc!=65535;
      std::string k=std::to_string(key);
      if(!counts.contains(k)) counts[k]={{"seen",0},{"clean",0},{"above_threshold",0}};
      counts[k]["seen"]=counts[k]["seen"].get<Long64_t>()+1;
      if(clean)counts[k]["clean"]=counts[k]["clean"].get<Long64_t>()+1;
      if(selected_ring(key)) {
        hist.at("ge"+k+"_all")->Fill(e);
        if(clean)hist.at("ge"+k+"_clean")->Fill(e);
      }
      if(e<=20)continue;
      counts[k]["above_threshold"]=counts[k]["above_threshold"].get<Long64_t>()+1;
      Long64_t time=failed?ts*8:Long64_t((ts*2-cfds+cfd/16384.)*4.);
      time=Long64_t(time+offsets.at(key));
      hits.push_back({time,i,e,key,clean});
      if(i%10000000==0)std::cout<<"read "<<i<<" / "<<tree->GetEntries()<<std::endl;
    }
    std::sort(hits.begin(),hits.end(),[](const auto& a,const auto& b) {
      return std::tie(a.time,a.key,a.entry)<std::tie(b.time,b.key,b.entry);
    });
    const auto before_dedup=hits.size();
    auto last=std::unique(hits.begin(),hits.end(),[](const auto& a,const auto& b) {
      return a.time==b.time&&a.key==b.key;
    });
    hits.erase(last,hits.end());
    json grouping=json::object();
    struct Mode{std::string name;int window;bool clean,backward;};
    for(const Mode& mode:std::vector<Mode>{{"legacy200",200,false,true},{"legacy100",100,false,true},
                                         {"clean200",200,true,true},{"forward200",200,false,false}}) {
      for(const std::string g:{"all","ring3","ring4","ring5"})make(mode.name+"_"+g);
      std::vector<size_t> kept;kept.reserve(hits.size());
      for(size_t i=0;i<hits.size();++i)if(!mode.clean||hits[i].clean)kept.push_back(i);
      Long64_t groups=0,members=0,pairs=0,repeat_pairs=0,rejected=0,emitted=0;
      size_t begin=0;
      for(size_t pos=0;pos<kept.size();) {
        const auto& ref=hits[kept[pos]];
        while(begin<pos&&hits[kept[begin]].time<ref.time-mode.window)++begin;
        size_t end=pos+1;
        while(end<kept.size()&&hits[kept[end]].time<=ref.time+mode.window)++end;
        const size_t first=mode.backward?begin:pos;
        ++groups;members+=end-first;
        std::vector<size_t> clover;
        auto emit=[&](int key,double energy) {
          int ring=selected_ring(key);if(!ring)return;
          hist.at(mode.name+"_all")->Fill(energy);
          hist.at(mode.name+"_ring"+std::to_string(ring))->Fill(energy);++emitted;
        };
        for(size_t i=first;i<end;++i) {
          const auto& h=hits[kept[i]];
          if(h.key>=212&&h.key<=215)clover.push_back(kept[i]);
          else emit(h.key,h.energy);
        }
        if(clover.size()==1)emit(hits[clover[0]].key,hits[clover[0]].energy);
        else if(clover.size()==2&&std::abs(hits[clover[0]].key-hits[clover[1]].key)!=2) {
          const auto& a=hits[clover[0]];const auto& b=hits[clover[1]];
          emit(a.energy>b.energy?a.key:b.key,a.energy+b.energy);++pairs;
          if(a.key==b.key)++repeat_pairs;
        } else if(!clover.empty())++rejected;
        pos=end;
      }
      grouping[mode.name]={{"references",groups},{"members",members},{"unique_input",kept.size()},
        {"selected_emitted",emitted},{"clover_pairs",pairs},{"same_crystal_pairs",repeat_pairs},
        {"clover_rejected_groups",rejected},{"window_ns",mode.window},{"backward_reuse",mode.backward}};
      std::cout<<mode.name<<" "<<grouping[mode.name].dump()<<std::endl;
    }
    require(original.same(Stat(input))&&text==read_text(conf),"Source/config changed");
    Digest digest;std::ifstream source(input,std::ios::binary);std::vector<char> block(8*1024*1024);
    while(source){source.read(block.data(),block.size());auto n=source.gcount();if(n>0)digest.add(block.data(),size_t(n));}
    require(source.eof()&&original.same(Stat(input)),"Source hashing failed");
    json meta={{"schema","sm150_legacy_efficiency_recovery_v1"},{"input",input.string()},
      {"input_sha256",digest.finish()},{"config_sha256",hash(text)},{"config",config},
      {"entries_read",tree->GetEntries()},{"channel_counts",counts},{"grouping",grouping},
      {"duplicate_time_channel_records_removed",before_dedup-hits.size()},
      {"elapsed_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()},
      {"exact_historical_reproduction",false},
      {"limits","Source energy file, source-only build variant, RNG seed and fit ranges were not preserved; reconstructed assumptions are explicit. Final singleton is retained. No particle gate, random subtraction or Doppler."}};
    fs::create_directories(output.parent_path());
    std::unique_ptr<TFile> out(TFile::Open(output.c_str(),"CREATE","Recovered source spectra",404));
    require(out&&!out->IsZombie(),"Cannot create output");
    for(const auto& kv:hist)require(kv.second->Write()>0,"Histogram write failure");
    require(TObjString(meta.dump().c_str()).Write("metadata")>0,"Metadata write failure");
    out->Close();require(!out->TestBit(TFile::kWriteError),"Output write failure");
    std::ofstream receipt(output.string()+".json");receipt<<meta.dump(2)<<'\n';
    std::cout<<"complete "<<meta["elapsed_seconds"]<<" s"<<std::endl;
  }catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';return 1;}
}
