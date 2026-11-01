// Source-run spectra and explicit Clover addback diagnostics.
// No in-beam gain, particle coincidence, random subtraction or Doppler correction.
#include "sm150_event_calibration.h"
#include <TFile.h>
#include <TH1D.h>
#include <TObjString.h>
#include <TTree.h>
#include <TLeaf.h>
#include <chrono>
#include <iostream>
#include <memory>

using namespace sm150::diagnostic;
struct Hit { double time, energy; int key; bool clean; };
int main(int argc, char** argv) {
  try {
    require(argc==4,"Usage: source_spectra INPUT_C1.root ENERGY.json NEW_OUTPUT.root");
    const fs::path input=fs::canonical(argv[1]), calibration=fs::canonical(argv[2]);
    const fs::path output=fs::absolute(argv[3]), receipt=output.string()+".json";
    require(!fs::exists(output) && !fs::exists(receipt),"Refusing existing output");
    const auto start=std::chrono::steady_clock::now();
    const Stat before(input);
    const std::string config_text=read_text(calibration);
    const json config=json::parse(config_text);
    require(config.at("run")=="0722","Expected source-run calibration");
    std::map<int,std::array<double,3>> coeff;
    for(auto it=config.at("channels").begin();it!=config.at("channels").end();++it)
      coeff.emplace(std::stoi(it.key()),it.value().at("coefficients").get<std::array<double,3>>());
    const std::map<int,std::vector<int>> groups={{3,{204,205,206,212,213,214,215}},
        {4,{211,304,305,306,308,507,508,510}},{5,{209,310,311,506,509}}};
    std::map<int,int> rings;
    for(const auto& g:groups)for(int key:g.second)rings[key]=g.first;
    require(coeff.size()==20,"Expected all 20 selected channels");
    for(const auto& kv:rings) require(coeff.count(kv.first),"Missing source coefficient");
    std::unique_ptr<TFile> fi(TFile::Open(input.c_str(),"READ"));
    require(fi && !fi->IsZombie() && !fi->TestBit(TFile::kRecovered),"Invalid source ROOT");
    auto* tree=dynamic_cast<TTree*>(fi->Get("tree")); require(tree,"Missing tree");
    tree->SetBranchStatus("*",false);
    const std::map<std::string,std::string> types={{"cid","Short_t"},{"sid","Short_t"},{"ch","Short_t"},
      {"evte","UShort_t"},{"ts","Long64_t"},{"sr","Short_t"},{"pileup","Bool_t"},{"outofr","Bool_t"}};
    for(const auto& kv:types) {
      auto* leaf=tree->GetLeaf(kv.first.c_str());
      require(leaf && kv.second==leaf->GetTypeName() && !leaf->GetLeafCount() && leaf->GetLenStatic()==1,"Invalid branch "+kv.first);
      tree->SetBranchStatus(kv.first.c_str(),true);
    }
    Short_t cid=0,sid=0,ch=0,sr=0;UShort_t adc=0;Long64_t ts=0;Bool_t pileup=false,outofr=false;
    tree->SetBranchAddress("cid",&cid);tree->SetBranchAddress("sid",&sid);tree->SetBranchAddress("ch",&ch);
    tree->SetBranchAddress("evte",&adc);tree->SetBranchAddress("ts",&ts);tree->SetBranchAddress("sr",&sr);
    tree->SetBranchAddress("pileup",&pileup);tree->SetBranchAddress("outofr",&outofr);
    tree->SetCacheSize(16*1024*1024);
    for(const auto& kv:types)tree->AddBranchToCache(kv.first.c_str(),true);
    tree->StopCacheLearningPhase();
    std::map<std::string,std::unique_ptr<TH1D>> hist;
    auto make=[&](const std::string& name,int bins=4096,double lo=0,double hi=4096) {
      auto p=std::make_unique<TH1D>(name.c_str(),";Energy (keV);Counts / bin",bins,lo,hi);
      p->SetDirectory(nullptr);hist[name]=std::move(p);
    };
    for(const auto& kv:coeff)make("ge"+std::to_string(kv.first));
    for(const std::string& group:{"all","ring3","ring4","ring5"})make(group+"_singles");
    for(int w:{50,100,200}) {make("clover2_ab"+std::to_string(w));make("ring3_ab"+std::to_string(w));make("all_ab"+std::to_string(w));}
    for(int a=212;a<=215;++a)for(int b=a+1;b<=215;++b)make("dt_"+std::to_string(a)+"_"+std::to_string(b),500,-2000,2000);
    std::vector<Hit> clover;clover.reserve(7000000);
    std::map<int,uint64_t> seen,clean,accepted;
    Long64_t previous=std::numeric_limits<Long64_t>::min();
    for(Long64_t i=0;i<tree->GetEntries();++i) {
      require(tree->GetEntry(i)>0,"Input read failure");
      require(sr==250 && cid==0,"Unexpected source module");
      require(ts>=previous,"Source coarse timestamp decreased");previous=ts;
      const int key=cid*10000+sid*100+ch;auto it=coeff.find(key);if(it==coeff.end())continue;
      ++seen[key];const bool good=!pileup && !outofr && adc!=65535;
      if(good)++clean[key];
      const auto& p=it->second;const double a=adc,e=p[0]+p[1]*a+p[2]*a*a;
      require(std::isfinite(e),"Nonfinite calibrated source energy");
      if(key>=212 && key<=215 && e>20)clover.push_back({double(ts)*8.,e,key,good});
      if(!good)continue;
      ++accepted[key];hist.at("ge"+std::to_string(key))->Fill(e);
      hist.at("all_singles")->Fill(e);hist.at("ring"+std::to_string(rings.at(key))+"_singles")->Fill(e);
      if(key<212||key>215)for(int w:{50,100,200}) {
        hist.at("all_ab"+std::to_string(w))->Fill(e);
        if(rings.at(key)==3)hist.at("ring3_ab"+std::to_string(w))->Fill(e);
      }
    }
    json addback=json::object();
    for(size_t i=0;i<clover.size();++i)for(size_t j=i+1;j<clover.size() && clover[j].time-clover[i].time<2000;++j) {
      const auto& a=clover[i];const auto& b=clover[j];
      if(a.key==b.key || !a.clean || !b.clean)continue;
      const int lo=std::min(a.key,b.key),hi=std::max(a.key,b.key);
      const double dt=a.key==lo?b.time-a.time:a.time-b.time;
      hist.at("dt_"+std::to_string(lo)+"_"+std::to_string(hi))->Fill(dt);
    }
    for(int window:{50,100,200}) {
      uint64_t one=0,two=0,rejected=0,bad_quality=0,same_channel=0;
      for(size_t i=0;i<clover.size();) {
        size_t j=i+1;while(j<clover.size() && clover[j].time-clover[i].time<=window)++j;
        bool good=true;double e=0;
        for(size_t k=i;k<j;++k) {good=good&&clover[k].clean;e+=clover[k].energy;}
        const bool single=j-i==1;
        const bool pair=j-i==2 && clover[i].key!=clover[i+1].key && std::abs(clover[i].key-clover[i+1].key)!=2;
        if(j-i==2 && clover[i].key==clover[i+1].key)++same_channel;
        if(!good)++bad_quality;
        else if(single||pair) {
          if(single)++one;else ++two;
          hist.at("clover2_ab"+std::to_string(window))->Fill(e);
          hist.at("ring3_ab"+std::to_string(window))->Fill(e);hist.at("all_ab"+std::to_string(window))->Fill(e);
        } else ++rejected;
        i=j;
      }
      addback[std::to_string(window)]={{"single_groups",one},{"adjacent_pairs",two},{"rejected_topology",rejected},
        {"rejected_quality",bad_quality},{"same_channel_pairs",same_channel}};
    }
    require(before.same(Stat(input)),"Source changed");
    require(read_text(calibration)==config_text,"Calibration changed");
    // Hash the source in bounded memory and retain the fresh checksum.
    Digest digest;std::ifstream source(input,std::ios::binary);std::vector<char> block(8*1024*1024);
    while(source) {source.read(block.data(),block.size());auto n=source.gcount();if(n>0)digest.add(block.data(),size_t(n));}
    require(source.eof() && before.same(Stat(input)),"Source hash failed or changed");
    json meta={{"schema","sm150_ge_source_spectra_v1"},{"input",input.string()},{"input_sha256",digest.finish()},
       {"entries_read",tree->GetEntries()},{"energy_config_sha256",hash(config_text)},{"run","0722"},
       {"seen",seen},{"clean",clean},{"accepted",accepted},{"addback",addback},
       {"source_timing","native coarse 8*ts; common module; no in-beam timing profile or CFD selection"},
       {"addback_rule","Clover2 only, forward cluster from earliest hit; >20keV per crystal; one crystal or two distinct adjacent crystals; other topologies rejected"},
       {"comparison_limit","Source diagnostic; not exact production grouping, which uses all Ge and in-beam corrected time"},
       {"quality","reject pileup OR outofr OR ADC=65535; normal and failed CFD combined"},
       {"elapsed_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()}};
    fs::create_directories(output.parent_path());
    std::unique_ptr<TFile> fo(TFile::Open(output.c_str(),"CREATE","Source efficiency diagnostics",404));
    require(fo && !fo->IsZombie(),"Cannot create output");
    for(const auto& item:hist)require(item.second->Write()>0,"Histogram write failed");
    require(TObjString(meta.dump().c_str()).Write("metadata")>0,"Metadata write failed");
    fo->Close();require(!fo->TestBit(TFile::kWriteError),"Output write failed");
    std::ofstream out(receipt);out<<meta.dump(2)<<'\n';out.close();require(bool(out),"Receipt write failed");
    std::cout<<meta.dump()<<std::endl;
  } catch(const std::exception& e) {std::cerr<<"source_spectra: "<<e.what()<<'\n';return 1;}
}
