// Keep prompt/random energy spectra separately for each ring-sector pixel.
// All event gates are independent of front/back energy agreement.
#include "sm150_event_calibration.h"
#include <TFile.h>
#include <TObjString.h>
#include <TH3I.h>
#include <TTree.h>
#include <array>
#include <chrono>
#include <iostream>
#include <memory>
#include <unistd.h>
using namespace sm150::diagnostic;
int main(int argc,char** argv) {
  try {
    fs::path input,output;
    for(int i=1;i<argc;++i) {
      std::string option=argv[i];
      auto value=[&](){require(++i<argc,"Missing argument");return std::string(argv[i]);};
      if(option=="--input")input=value();else if(option=="--output")output=value();
      else throw std::runtime_error("Unknown argument: "+option);
    }
    require(!input.empty()&&!output.empty()&&!fs::exists(output)&&!fs::exists(output.string()+".partial"),"Missing paths or existing output");
    const auto start=std::chrono::steady_clock::now();const Stat before(input);
    std::unique_ptr<TFile> source(TFile::Open(input.c_str(),"READ"));
    require(source&&!source->IsZombie(),"Cannot open input");
    auto* tree=dynamic_cast<TTree*>(source->Get("pairs"));
    auto* metadata=dynamic_cast<TObjString*>(source->Get("metadata"));require(tree&&metadata,"Missing input tree/metadata");
    json meta=json::parse(metadata->GetString().Data());require(meta.at("schema")=="sm150_s3_joint_sample_v1","Wrong input schema");
    // Historical samples have forty segments; new independent checks may use
    // more. Validate against the source receipt instead of hard-coding forty.
    const int blocks=meta.contains("segments")?int(meta.at("segments").size()):40;
    require(blocks>=4&&blocks<=1000,"Invalid source segment count");
    Int_t block=0,r=0,s=0,nr=0,ns=0,flags=0;Double_t re=0,se=0,dt=0;
    tree->SetBranchStatus("*",false);
    for(const auto* name:{"block","r","s","nr","ns","flags","r_e0_keV","s_e0_keV","dt_ns"})tree->SetBranchStatus(name,true);
    tree->SetBranchAddress("block",&block);tree->SetBranchAddress("r",&r);tree->SetBranchAddress("s",&s);
    tree->SetBranchAddress("nr",&nr);tree->SetBranchAddress("ns",&ns);tree->SetBranchAddress("flags",&flags);
    tree->SetBranchAddress("r_e0_keV",&re);tree->SetBranchAddress("s_e0_keV",&se);tree->SetBranchAddress("dt_ns",&dt);
    std::array<std::unique_ptr<TH3I>,4> spectra;
    const std::array<std::string,4> names={"ring_prompt","ring_random","sector_prompt","sector_random"};
    for(int i=0;i<4;++i) {
      spectra[i]=std::make_unique<TH3I>(names[i].c_str(),"E0;E0 [keV];pixel (ring*32+sector);block modulo4",800,0,200000,768,0,768,4,0,4);
      spectra[i]->SetDirectory(nullptr);
    }
    std::array<uint64_t,4> prompt{},random{};uint64_t clean=0;
    for(Long64_t i=0;i<tree->GetEntries();++i) {
      require(tree->GetEntry(i)>0,"Pair read failed");require(r>=0&&r<24&&s>=24&&s<56&&block>=0&&block<blocks,"Invalid pair labels");
      if(flags||nr!=1||ns!=1)continue;
      require(std::isfinite(re)&&std::isfinite(se)&&std::isfinite(dt),"Non-finite pair");++clean;
      int region=-1;
      if(dt>=-80&&dt<80)region=0;else if((dt>=-880&&dt< -720)||(dt>=720&&dt<880))region=1;
      if(region<0)continue;
      const int pixel=r*32+s-24,split=block%4;
      spectra[region]->Fill(re,pixel+.5,split+.5);spectra[2+region]->Fill(se,pixel+.5,split+.5);
      (region?random:prompt)[split]++;
    }
    require(before.same(Stat(input)),"Input changed");
    fs::create_directories(fs::absolute(output).parent_path());
    const fs::path partial=output.string()+".partial";
    std::unique_ptr<TFile> dest(TFile::Open(partial.c_str(),"CREATE","Conditional paired S3 spectra",404));require(dest&&!dest->IsZombie(),"Cannot create output");
    for(auto& h:spectra) {
      double total=h->Integral(0,801,0,769,0,5);
      const auto& counts=std::string(h->GetName()).find("prompt")!=std::string::npos?prompt:random;
      require(total==double(std::accumulate(counts.begin(),counts.end(),uint64_t(0))),"Histogram count loss");h->Write();
    }
    json receipt={{"schema","sm150_s3_conditioned_spectra_v1"},{"status","complete"},{"run",meta.at("run")},{"period",meta.at("period")},
      {"input",fs::absolute(input).string()},{"input_sha256",hash(read_text(input))},{"input_metadata",meta},
      {"selection","flags=0,nr=ns=1; no added energy or energy-difference gate; source thresholds and guards retained"},
      {"prompt_ns",{-80,80}},{"random_ns",{{-880,-720},{720,880}}},{"alpha",.5},{"clean_pairs",clean},
      {"prompt_by_split",prompt},{"random_by_split",random},{"bin_keV",250},{"maximum_energy_keV",200000},{"source_block_count",blocks},
      {"program_sha256",hash(read_text(fs::canonical("/proc/self/exe")))},
      {"elapsed_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()}};
    TObjString record(receipt.dump().c_str());record.Write("metadata");dest->Close();
    require(::link(partial.c_str(),output.c_str())==0&&::unlink(partial.c_str())==0,"Cannot publish output without overwriting");
    std::ofstream out(output.string()+".json");out<<receipt.dump(2)<<'\n';out.close();require(bool(out),"Cannot write receipt");
    std::cout<<meta.at("run")<<" complete clean="<<clean<<" output="<<output<<std::endl;return 0;
  }catch(const std::exception& e){std::cerr<<"conditioned S3 extraction: "<<e.what()<<'\n';return 1;}
}
