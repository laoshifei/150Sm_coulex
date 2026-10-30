// Extract distributed S3 singles; exclude attenuation and front/back corrections.
#include "sm150_event_calibration.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLeaf.h>
#include <TNamed.h>
#include <TTree.h>
#include <TUUID.h>
#include <chrono>
#include <iostream>
#include <memory>
#include <unistd.h>
using namespace sm150::diagnostic;

int main(int argc,char** argv) {
  try {
    int run=0,blocks=40;Long64_t sample=2000000;
    fs::path project=fs::canonical("/proc/self/exe").parent_path().parent_path(),input,output;
    for(int i=1;i<argc;++i) {
      std::string a=argv[i];
      if(a=="--help") {
        std::cout<<"extract RUN --output FILE [--input C2.root] [--project DIR]\n"
          "  --sample-entries N  Distributed sample; 0 reads all entries\n"
          "  --segments N       Equal source-entry intervals; default 40\n";return 0;
      }
      auto value=[&](){require(++i<argc,"Missing option value");return std::string(argv[i]);};
      auto number=[](const std::string& s){size_t pos=0;auto x=std::stoll(s,&pos);require(pos==s.size(),"Invalid integer");return x;};
      if(a=="--input")input=value();else if(a=="--output")output=value();
      else if(a=="--project")project=value();else if(a=="--sample-entries")sample=number(value());
      else if(a=="--segments")blocks=number(value());
      else if(!a.empty() && a[0]!='-' && !run)run=number(a);
      else throw std::runtime_error("Unknown argument: "+a);
    }
    require(run>0 && run<10000 && blocks>=4 && blocks<=200 && blocks%2==0 && sample>=0 && !output.empty(),"Invalid options; use --help");
    project=fs::canonical(project);
    if(input.empty())input="/mnt/e/Sm150/rootfile/data_C2_"+rn(run)+"_"+rn(run)+".root";
    input=fs::canonical(input);output=fs::absolute(output);
    fs::path partial=output.string()+".partial",receipt=output.string()+".json";
    require(!fs::exists(output) && !fs::exists(partial) && !fs::exists(receipt),"Output or partial exists");
    Parameters pars(project,run,false);pars.verify();const Stat before(input);
    const auto started=std::chrono::steady_clock::now();
    std::unique_ptr<TFile> source(TFile::Open(input.c_str(),"READ"));
    require(source && !source->IsZombie(),"Cannot open C2 input");
    auto* tree=dynamic_cast<TTree*>(source->Get("tree"));require(tree,"Missing decoded tree");
    const Long64_t total=tree->GetEntries(),wanted=sample?std::min(sample,total):total;
    require(total>0 && wanted>=blocks,"Too few entries for requested segments");
    Short_t cid=0,sid=0,ch=0;UShort_t adc=0;Bool_t pileup=false,outofr=false;
    const std::map<std::string,std::string> types={{"cid","Short_t"},{"sid","Short_t"},{"ch","Short_t"},{"evte","UShort_t"},{"pileup","Bool_t"},{"outofr","Bool_t"}};
    tree->SetBranchStatus("*",false);
    for(const auto& kv:types) {
      auto* leaf=tree->GetLeaf(kv.first.c_str());require(leaf && leaf->GetTypeName()==kv.second,"Wrong decoded branch: "+kv.first);
      tree->SetBranchStatus(kv.first.c_str(),true);
    }
    require(tree->SetBranchAddress("cid",&cid)>=0 && tree->SetBranchAddress("sid",&sid)>=0 &&
      tree->SetBranchAddress("ch",&ch)>=0 && tree->SetBranchAddress("evte",&adc)>=0 &&
      tree->SetBranchAddress("pileup",&pileup)>=0 && tree->SetBranchAddress("outofr",&outofr)>=0,"Cannot bind decoded branches");
    tree->SetCacheSize(8*1024*1024);
    for(const auto& kv:types)tree->AddBranchToCache(kv.first.c_str(),true);
    tree->StopCacheLearningPhase();
    fs::create_directories(output.parent_path());
    std::unique_ptr<TFile> target(TFile::Open(partial.c_str(),"CREATE","S3 energy before attenuation",404));
    require(target && !target->IsZombie(),"Cannot create isolated output");target->cd();
    std::array<TH2D*,56> clean{};std::array<TH1D*,56> inclusive{};
    json channel_records=json::array();std::map<int,int> indices;
    for(const auto& kv:pars.channels) {
      const auto& c=kv.second;if(c.kind<3)continue;
      const int index=c.kind==3?c.ring-1:24+c.sector-1;
      const auto name=(c.kind==3?"ring":"sector")+rn(c.kind==3?c.ring:c.sector).substr(2);
      indices[kv.first]=index;
      clean[index]=new TH2D(("clean_"+name).c_str(),"Clean S3 singles;Energy before attenuation [keV];Source segment",3000,0,150000,blocks,0,blocks);
      inclusive[index]=new TH1D(("inclusive_"+name).c_str(),"All S3 singles;Energy before attenuation [keV];Counts",3000,0,150000);
      channel_records.push_back({{"index",index},{"name",name},{"key",kv.first},{"old_k",c.att},
        {"calibration",c.cal},{"hardware_gain",c.gain},{"k_ge",c.kg}});
    }
    require(indices.size()==56,"Incomplete S3 channel map");
    json segments=json::array();Long64_t processed=0,s3_count=0,clean_count=0;
    std::array<Long64_t,56> hits{},accepted{};
    for(int block=0;block<blocks;++block) {
      const Long64_t low=total*block/blocks,high=total*(block+1)/blocks;
      const Long64_t count=wanted*(block+1)/blocks-wanted*block/blocks;
      const Long64_t first=low+(high-low-count)/2,last=first+count;
      require(first>=low && last<=high,"Invalid sample interval");
      tree->SetCacheEntryRange(first,last);Long64_t ns3=0,nclean=0;
      for(Long64_t entry=first;entry<last;++entry) {
        require(tree->GetEntry(entry)>0,"Failed event read");++processed;
        if(cid!=1)continue;
        const int key=10000+100*sid+ch;auto found=indices.find(key);if(found==indices.end())continue;
        const int index=found->second;const auto& c=pars.channels.at(key);
        const double x=adc,e=(c.cal[0]+c.cal[1]*x+c.cal[2]*x*x)*c.gain*c.kg;
        require(std::isfinite(e),"Nonfinite energy");
        ++s3_count;++ns3;++hits[index];inclusive[index]->Fill(e);
        if(pileup || outofr || adc==65535)continue;
        clean[index]->Fill(e,block+0.5);++clean_count;++nclean;++accepted[index];
      }
      segments.push_back({{"block",block},{"source_interval",{low,high}},{"read_interval",{first,last}},
        {"entries",count},{"s3_hits",ns3},{"clean_hits",nclean},{"split",block%2?"validation":"fit"}});
      const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
      std::cout<<"progress run="<<rn(run)<<" segment="<<block+1<<'/'<<blocks<<" read="<<processed<<" clean="<<clean_count<<" seconds="<<seconds<<std::endl;
    }
    require(processed==wanted,"Incomplete extraction");pars.verify();require(before.same(Stat(input)),"Input changed");
    target->cd();
    for(int i=0;i<56;++i) {
      require(Long64_t(clean[i]->GetEntries())==accepted[i] && Long64_t(inclusive[i]->GetEntries())==hits[i],"Histogram count mismatch");
      clean[i]->Write();inclusive[i]->Write();
    }
    json hashes=json::object();for(auto it=pars.snapshots.begin();it!=pars.snapshots.end();++it)hashes[it.key()]=it.value().at("sha256");
    json metadata={{"schema","sm150_s3_attenuation_spectra_v1"},{"status","complete"},{"run",rn(run)},
      {"period",pars.group},{"input",input.string()},{"input_bytes",before.value.st_size},{"input_uuid",source->GetUUID().AsString()},
      {"input_entries",total},{"sample_entries",processed},{"sample_fraction",double(processed)/total},
      {"full_input_read",processed==total},{"segments",segments},{"channels",channel_records},{"parameter_hashes",hashes},
      {"s3_hits",s3_count},{"clean_hits",clean_count},{"channel_hits",hits},{"channel_clean_hits",accepted},
      {"energy_definition","ADC polynomial times hardware gain; no attenuation or front/back normalization"},
      {"selection","S3 singles; clean excludes pileup/outofr/ADC65535; no time, Ge, multiplicity or front/back energy gate"},
      {"energy_bins",3000},{"energy_low_keV",0},{"energy_high_keV",150000},
      {"elapsed_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()},
      {"program_sha256",hash(read_text(fs::canonical("/proc/self/exe")))}};
    TNamed provenance("metadata",metadata.dump().c_str());provenance.Write();target->Close();source->Close();
    std::unique_ptr<TFile> check(TFile::Open(partial.c_str(),"READ"));
    require(check && !check->IsZombie() && check->Get("metadata") && check->Get("clean_sector32"),"Cannot reopen output");check->Close();
    require(::link(partial.c_str(),output.c_str())==0,"Cannot publish without overwriting");
    require(::unlink(partial.c_str())==0,"Cannot remove partial alias");
    std::ofstream report(receipt);report<<metadata.dump(2)<<'\n';report.close();require(bool(report),"Cannot save receipt");
    std::cout<<"complete run="<<rn(run)<<" read="<<processed<<" clean="<<clean_count<<" output="<<output<<std::endl;return 0;
  }catch(const std::exception& e){std::cerr<<"s3 attenuation extraction: "<<e.what()<<'\n';return 1;}
}
