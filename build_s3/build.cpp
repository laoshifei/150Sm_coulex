// S3 calibration: decoded ADC -> frozen timing -> unmerged temporal clusters.
// No Ge gate or ring-sector energy-agreement cut is applied to this source.
#include "build.h"
#include "sm150_event_calibration.h"
#include "sm150_time_alignment.h"
#include <TFile.h>
#include <TLeaf.h>
#include <TObjString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <limits>
#include <memory>
#include <queue>
#include <unistd.h>

namespace {
using namespace sm150::diagnostic;
using Clock=std::chrono::steady_clock;
struct Hit {
  Long64_t source=0,native=0,time=0;
  int key=0,face=0,id=0,adc=0,flags=0;
  double energy=0.,previous=0.;
};
struct Later {
  bool operator()(const Hit& a,const Hit& b) const {
    return a.time>b.time || (a.time==b.time && a.source>b.source);
  }
};
struct Output {
  TTree events{"s3_events","Unmerged S3 clusters; no energy-agreement gate"};
  TTree pairs{"s3_pairs","Exactly one raw ring and sector in a temporal cluster"};
  Long64_t cluster=0,first_entry=0,ring_entry=0,sector_entry=0;
  Long64_t ring_native=0,sector_native=0,ring_time=0,sector_time=0;
  Int_t block=0,n_ring=0,n_sector=0,ring=0,sector=0,ring_key=0,sector_key=0;
  Int_t ring_adc=0,sector_adc=0,quality=0;
  Bool_t boundary=false;
  Double_t ring_energy=0,sector_energy=0,ring_previous=0,sector_previous=0,dt=0;
  std::vector<int> keys,faces,ids,adcs,flags;
  std::vector<Long64_t> entries,natives,times;
  std::vector<double> energies,previous;
  uint64_t clusters=0,paired_clusters=0,single_face_clusters=0,pair_count=0,good_pairs=0,emitted_hits=0;
  size_t max_multiplicity=0;
  Output() {
    events.Branch("cluster",&cluster);events.Branch("block",&block);
    events.Branch("first_entry",&first_entry);events.Branch("boundary",&boundary);
    events.Branch("n_ring_raw",&n_ring);events.Branch("n_sector_raw",&n_sector);
    events.Branch("channel",&keys);events.Branch("face",&faces);events.Branch("id",&ids);
    events.Branch("adc",&adcs);events.Branch("flags",&flags);events.Branch("source_entry",&entries);
    events.Branch("native_time_ns",&natives);events.Branch("time_ns",&times);
    events.Branch("energy_pre_fb_keV",&energies);events.Branch("energy_previous_keV",&previous);
    pairs.Branch("cluster",&cluster);pairs.Branch("block",&block);pairs.Branch("boundary",&boundary);
    pairs.Branch("ring_id",&ring);pairs.Branch("sector_id",&sector);
    pairs.Branch("ring_key",&ring_key);pairs.Branch("sector_key",&sector_key);
    pairs.Branch("ring_adc",&ring_adc);pairs.Branch("sector_adc",&sector_adc);
    pairs.Branch("ring_entry",&ring_entry);pairs.Branch("sector_entry",&sector_entry);
    pairs.Branch("ring_native_ns",&ring_native);pairs.Branch("sector_native_ns",&sector_native);
    pairs.Branch("ring_time_ns",&ring_time);pairs.Branch("sector_time_ns",&sector_time);
    pairs.Branch("quality_flags",&quality);pairs.Branch("dt_sector_minus_ring_ns",&dt);
    pairs.Branch("ring_energy_pre_fb_keV",&ring_energy);pairs.Branch("sector_energy_pre_fb_keV",&sector_energy);
    pairs.Branch("ring_energy_previous_keV",&ring_previous);pairs.Branch("sector_energy_previous_keV",&sector_previous);
    events.SetAutoFlush(-16*1024*1024);pairs.SetAutoFlush(-16*1024*1024);
  }
  void flush(const std::vector<Hit>& hits,bool edge,int64_t limit,int blocks) {
    if(hits.empty())return;
    cluster=clusters++;boundary=edge;first_entry=hits.front().source;
    n_ring=n_sector=0;keys.clear();faces.clear();ids.clear();adcs.clear();flags.clear();
    entries.clear();natives.clear();times.clear();energies.clear();previous.clear();
    max_multiplicity=std::max(max_multiplicity,hits.size());
    const Hit *r=nullptr,*s=nullptr;
    for(const auto& h:hits) {
      ++emitted_hits;first_entry=std::min(first_entry,h.source);
      if(h.face==0){++n_ring;r=&h;}else{++n_sector;s=&h;}
      keys.push_back(h.key);faces.push_back(h.face);ids.push_back(h.id);adcs.push_back(h.adc);
      flags.push_back(h.flags);entries.push_back(h.source);natives.push_back(h.native);times.push_back(h.time);
      energies.push_back(h.energy);previous.push_back(h.previous);
    }
    block=std::min(blocks-1,int((static_cast<long double>(first_entry)*blocks)/limit));
    if(!n_ring || !n_sector){++single_face_clusters;return;}
    ++paired_clusters;require(events.Fill()>0,"S3 cluster write failed");
    if(n_ring!=1 || n_sector!=1)return;
    ring=r->id;sector=s->id;ring_key=r->key;sector_key=s->key;
    ring_adc=r->adc;sector_adc=s->adc;ring_entry=r->source;sector_entry=s->source;
    ring_native=r->native;sector_native=s->native;ring_time=r->time;sector_time=s->time;
    ring_energy=r->energy;sector_energy=s->energy;ring_previous=r->previous;sector_previous=s->previous;
    quality=r->flags|(s->flags<<4);dt=double(s->time-r->time);
    require(pairs.Fill()>0,"S3 pair write failed");++pair_count;if(!quality && !boundary)++good_pairs;
  }
};
}

int build::Process() {
  using namespace sm150::diagnostic;
  const auto started=Clock::now();
  try {
    const auto& o=options_;
    fs::path project=o.project.empty()?fs::canonical("/proc/self/exe").parent_path().parent_path():fs::canonical(o.project);
    const auto input=fs::canonical(o.input.empty()?fs::path("/mnt/e/Sm150/rootfile/data_C2_"+rn(o.run)+"_"+rn(o.run)+".root"):o.input);
    const auto output=fs::absolute(o.output.empty()?project/"rootfile"/("data"+rn(o.run)+"_s3_calibration.root"):o.output);
    const fs::path partial=output.string()+".partial";
    require(!fs::exists(output) && !fs::exists(partial),"Output or unfinished output exists: "+output.string());
    require(!fs::exists(output.string()+".json"),"Existing extraction receipt");
    Parameters pars(project,o.run,true);
    fs::path profile=o.profile_dir;
    if(profile.empty() && pars.paths.active()) {
      profile=pars.paths.timing_directory();
    } else if(profile.empty()) {
      const auto active=project/"pars/time_alignment/active.txt";pars.freeze(active);
      std::istringstream in(read_text(active));std::string name,extra;
      require(bool(in>>name) && !(in>>extra) && name.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_")==std::string::npos,"Invalid active timing profile");
      profile=project/"pars/time_alignment/versions"/name;
    }
    profile=fs::canonical(profile);
    for(const auto& entry:fs::directory_iterator(profile))if(entry.is_regular_file() && entry.path().filename()!="README.md")pars.freeze(entry.path());
    sm150::TimeAlignment timing;timing.LoadProfile(profile.string(),o.run);
    require(timing.total(),"S3 calibration requires a complete timing profile");
    pars.verify();const Stat input_stat(input);
    std::unique_ptr<TFile> source(TFile::Open(input.c_str(),"READ"));
    require(source && !source->IsZombie(),"Cannot read decoded ROOT input");
    auto* tree=dynamic_cast<TTree*>(source->Get("tree"));require(tree,"Missing decoded tree");
    const Long64_t total=tree->GetEntries(),limit=o.entries?std::min<Long64_t>(o.entries,total):total;
    require(limit>0,"Empty input range");
    const std::map<std::string,std::string> types={{"cid","Short_t"},{"sid","Short_t"},{"ch","Short_t"},{"sr","Short_t"},{"evte","UShort_t"},{"ts","Long64_t"},{"pileup","Bool_t"},{"outofr","Bool_t"}};
    for(const auto& kv:types) {
      auto* leaf=tree->GetLeaf(kv.first.c_str());require(leaf && kv.second==leaf->GetTypeName(),"Wrong decoded branch: "+kv.first);
    }
    tree->SetCacheSize(32*1024*1024);
    for(const auto& kv:types)tree->AddBranchToCache(kv.first.c_str(),true);
    tree->StopCacheLearningPhase();
    fs::create_directories(output.parent_path());
    std::unique_ptr<TFile> target(TFile::Open(partial.c_str(),"CREATE","S3 normalization calibration source",404));
    require(target && !target->IsZombie(),"Cannot create independent calibration output");target->cd();
    Output out;
    TTreeReader reader(tree);
    TTreeReaderValue<Short_t> cid(reader,"cid"),sid(reader,"sid"),ch(reader,"ch"),sr(reader,"sr");
    TTreeReaderValue<UShort_t> adc(reader,"evte");TTreeReaderValue<Long64_t> timestamp(reader,"ts");
    TTreeReaderValue<Bool_t> pileup(reader,"pileup"),outofr(reader,"outofr");
    reader.SetEntriesRange(0,limit);
    std::priority_queue<Hit,std::vector<Hit>,Later> buffer;std::vector<Hit> cluster;
    Long64_t max_native=std::numeric_limits<Long64_t>::min(),last_time=max_native;
    uint64_t processed=0,s3_hits=0,accepted=0,energy_rejected=0,flagged=0;
    size_t peak_buffer=0;
    auto consume=[&](const Hit& h) {
      require(h.time>=last_time,"Hit exceeds timestamp reorder buffer");last_time=h.time;
      if(!cluster.empty() && h.time-cluster.back().time>o.cluster_gap_ns) {
        out.flush(cluster,out.clusters==0,limit,o.blocks);cluster.clear();
      }
      cluster.push_back(h);require(cluster.size()<1000000,"Unbounded S3 temporal cluster");
    };
    auto drain=[&](Long64_t watermark) {
      while(!buffer.empty() && buffer.top().time<=watermark) {auto h=buffer.top();buffer.pop();consume(h);}
    };
    while(reader.Next()) {
      ++processed;
      if(processed%5000000==0) {
        double seconds=std::chrono::duration<double>(Clock::now()-started).count();
        std::cout<<"progress run="<<rn(o.run)<<" entries="<<processed<<"/"<<limit<<" pairs="<<out.pair_count<<" seconds="<<seconds<<std::endl;
      }
      if(*cid!=1)continue;
      const int slot=*sid;
      if(slot!=11 && slot!=12 && slot!=(pars.group=="low"?2:9) && slot!=(pars.group=="low"?3:10))continue;
      const int key=10000+slot*100+*ch;auto found=pars.channels.find(key);
      if(found==pars.channels.end() || found->second.kind<3)continue;
      ++s3_hits;const auto& c=found->second;const double energy=c.energy(double(*adc),false);
      require(std::isfinite(energy),"Nonfinite S3 energy");
      if(energy<=o.min_energy_keV){++energy_rejected;continue;}
      const int step=*sr==250?8:*sr==100?10:0;require(step && *timestamp>=0 && *timestamp<std::numeric_limits<Long64_t>::max()/step,"Invalid decoded time");
      const auto native=*timestamp*step;
      const double delta=timing.TotalParticleOffsetFromAdc(key,*adc);
      require(std::isfinite(delta) && std::abs(delta)<o.reorder_ns/2.,"Timing correction exceeds reorder assumptions");
      if(max_native!=std::numeric_limits<Long64_t>::min())require(native>=max_native-o.reorder_ns/2,"Native S3 stream is not locally ordered");
      max_native=std::max(max_native,native);
      Hit h;h.source=reader.GetCurrentEntry();h.native=native;
      h.time=static_cast<Long64_t>(double(native)+delta);h.key=key;h.face=c.kind-3;
      h.id=h.face?c.sector:c.ring;h.adc=*adc;h.energy=energy;h.previous=c.energy(double(*adc),true);
      h.flags=(*pileup?1:0)|(*outofr?2:0)|(*adc==65535?4:0);if(h.flags)++flagged;
      buffer.push(h);++accepted;peak_buffer=std::max(peak_buffer,buffer.size());
      drain(max_native-o.reorder_ns);
    }
    require(processed==uint64_t(limit),"Input ended before requested range");
    drain(std::numeric_limits<Long64_t>::max());out.flush(cluster,true,limit,o.blocks);
    require(out.emitted_hits==accepted && out.clusters==out.paired_clusters+out.single_face_clusters,"S3 count conservation failed");
    pars.verify();require(input_stat.same(Stat(input)),"Input changed while reading");
    json parameter_hashes=json::object();
    for(auto it=pars.snapshots.begin();it!=pars.snapshots.end();++it)parameter_hashes[it.key()]=it.value().at("sha256");
    json metadata={{"schema","sm150_s3_calibration_v1"},{"status","complete"},{"run",rn(o.run)},{"period",pars.group},
      {"input",input.string()},{"input_uuid",source->GetUUID().AsString()},{"input_bytes",input_stat.value.st_size},
      {"input_entries",total},{"processed_entries",processed},{"full_input_read",limit==total},{"output",output.string()},
      {"timing_profile",timing.profile()},{"timing_source_sha256",timing.source_sha256()},
      {"timing_profile_dir",profile.string()},{"parameter_hashes",parameter_hashes},{"channels",pars.records},
      {"energy_definition","keV after calibration/gain/attenuation; pre_fb excludes previous normalization"},
      {"timestamp_definition","Native sr clock plus complete frozen ADC-dependent correction once; no dither"},
      {"cluster_definition","Connected components in corrected time; no merging/addback"},
      {"cluster_gap_ns",o.cluster_gap_ns},{"reorder_ns",o.reorder_ns},{"min_pre_fb_energy_keV",o.min_energy_keV},
      {"energy_agreement_gate",false},{"addback",false},{"blocks",o.blocks},
      {"block_definition","Equal source-entry blocks; even fit, odd validation; file-edge clusters flagged"},
      {"quality_flags","bits0/1/2 ring pileup/outofr/ADC saturation; sector shifted by4"},
      {"s3_hits",s3_hits},{"accepted_s3_hits",accepted},{"energy_rejected",energy_rejected},{"flagged_hits",flagged},
      {"clusters",out.clusters},{"clusters_both_faces",out.paired_clusters},{"clusters_single_face",out.single_face_clusters},
      {"single_pairs",out.pair_count},{"good_nonboundary_pairs",out.good_pairs},{"max_multiplicity",out.max_multiplicity},
      {"peak_reorder_buffer",peak_buffer},{"elapsed_seconds",std::chrono::duration<double>(Clock::now()-started).count()}};
    target->cd();require(out.events.Write()>0 && out.pairs.Write()>0,"Output tree write failed");
    const auto text=metadata.dump();require(TObjString(text.c_str()).Write("metadata")>0,"Metadata write failed");
    out.events.SetDirectory(nullptr);out.pairs.SetDirectory(nullptr);target->Close();
    std::unique_ptr<TFile> check(TFile::Open(partial.c_str(),"READ"));require(check && !check->IsZombie(),"Output readback failed");
    auto* ct=dynamic_cast<TTree*>(check->Get("s3_pairs"));auto* et=dynamic_cast<TTree*>(check->Get("s3_events"));
    require(ct && et && ct->GetEntries()==Long64_t(out.pair_count) && et->GetEntries()==Long64_t(out.paired_clusters),"Readback counts differ");check->Close();
    require(::link(partial.c_str(),output.c_str())==0,"Cannot publish without overwriting output");
    require(::unlink(partial.c_str())==0,"Cannot remove published temporary link");
    std::ofstream receipt(output.string()+".json");receipt<<metadata.dump(2)<<'\n';receipt.close();require(bool(receipt),"Cannot write receipt");
    std::cout<<"complete run="<<rn(o.run)<<" entries="<<processed<<" pairs="<<out.pair_count<<" output="<<output<<std::endl;
    return 0;
  }catch(const std::exception& e){std::cerr<<"build_s3: "<<e.what()<<'\n';return 1;}
}
