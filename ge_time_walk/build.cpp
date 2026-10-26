// 顺序：冻结参数 -> 刻度/原生时间 -> 全候选配对 -> 双面检查 -> 分CFD时间谱。
#include "build.h"
#include "raw.h"
#include "sm150_event_calibration.h"
#include "sm150_time_alignment.h"
#include <TFile.h>
#include <TH2D.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TTree.h>
#include <chrono>
#include <deque>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <tuple>
using namespace sm150::diagnostic;
namespace {
const std::set<int> GE={204,205,206,209,211,212,213,214,215,304,305,306,308,310,311,506,507,508,509,510};
constexpr double SHIFT_BOUND=1000000.; // 已校准时间与原生coarse时间的严格运行时上界/ns。
#pragma pack(push,1)
struct InputRecord {
  int16_t sr;bool pileup,outofr;int16_t cid,sid,ch;uint16_t evte;int64_t ts;int16_t cfd;bool cfdft;int16_t cfds;
};
struct ParticleIdentityRecord {int64_t ts;uint16_t evte;int16_t sr,cid,sid,ch;uint8_t pileup,outofr;};
#pragma pack(pop)
static_assert(sizeof(InputRecord)==25 && sizeof(ParticleIdentityRecord)==20,"Unexpected raw hash layout");
struct Hit {
  int key=0,kind=0,ring=0,sector=0;
  int64_t entry=0,coarse=0,native=0,ts=0;
  double e=0,t=0;
  uint16_t adc=0;
  int16_t sr=0,cfd=0,cfds=0;
  bool failed=false,done=false;
};
struct PairRecord {
  int ge_key=0,particle_key=0,particle_kind=0,particle_ring=0,particle_sector=0,half=0;
  int64_t ge_entry=0,particle_entry=0,sector_entry=-1,ge_ts=0,particle_ts=0;
  int64_t ge_native_ns=0,particle_native_ns=0;
  double ge_energy_keV=0,particle_energy_keV=0,dt_current_ns=0,dt_native_ge_ns=0,dt_native_both_ns=0;
  double particle_corrected_ns=0,ge_corrected_ns=0;
  uint16_t ge_adc=0,particle_adc=0;
  int16_t ge_cfd=0,ge_cfds=0,particle_cfd=0,particle_cfds=0,ge_sr=0,particle_sr=0;
  bool ge_failed=false,particle_failed=false,s3_unique_match=false;
  bool ge_pileup=false,ge_outofr=false,particle_pileup=false,particle_outofr=false;
};
struct Writer {
  const WalkOptions& o;TFile& file;uint64_t split;
  std::map<std::string,std::unique_ptr<TH2D>> hist;
  std::map<std::string,uint64_t> counts;
  json descriptions=json::object();
  TTree* tree=nullptr;PairRecord row;
  Writer(const WalkOptions& opt,TFile& f,uint64_t s):o(opt),file(f),split(s) {
    if(!o.write_pairs) return;
    tree=new TTree("pairs","All selected raw-entry Ge-particle pairs; shared hits are correlated");
    tree->SetAutoFlush(-8*1024*1024);
#define BRANCH(x) tree->Branch(#x,&row.x)
    BRANCH(ge_key);BRANCH(particle_key);BRANCH(particle_kind);BRANCH(particle_ring);BRANCH(particle_sector);BRANCH(half);
    BRANCH(ge_entry);BRANCH(particle_entry);BRANCH(sector_entry);BRANCH(ge_ts);BRANCH(particle_ts);
    BRANCH(ge_native_ns);BRANCH(particle_native_ns);BRANCH(ge_energy_keV);BRANCH(particle_energy_keV);
    BRANCH(dt_current_ns);BRANCH(dt_native_ge_ns);BRANCH(dt_native_both_ns);BRANCH(particle_corrected_ns);BRANCH(ge_corrected_ns);
    BRANCH(ge_adc);BRANCH(particle_adc);BRANCH(ge_cfd);BRANCH(ge_cfds);BRANCH(particle_cfd);BRANCH(particle_cfds);
    BRANCH(ge_sr);BRANCH(particle_sr);BRANCH(ge_failed);BRANCH(particle_failed);BRANCH(s3_unique_match);
    BRANCH(ge_pileup);BRANCH(ge_outofr);BRANCH(particle_pileup);BRANCH(particle_outofr);
#undef BRANCH
  }
  TH2D& book(int ge,const std::string& ref,int state,int half,const std::string& coordinate) {
    std::string name="ge"+std::to_string(ge)+"_"+ref+(state?"_failed":"_normal")+"_half"+std::to_string(half)+"_"+coordinate;
    auto it=hist.find(name);if(it!=hist.end()) return *it->second;
    int ne=int(std::llround((o.ge_max-o.ge_min)/o.energy_bin_keV));
    int nt=int(std::llround(2*o.time_range_ns/o.time_bin_ns));
    auto h=std::make_unique<TH2D>(name.c_str(),(name+";Ge laboratory energy (keV);t_{particle}-t_{Ge} (ns)").c_str(),ne,o.ge_min,o.ge_max,nt,-o.time_range_ns,o.time_range_ns);
    h->SetDirectory(nullptr);
    descriptions[name]={{"ge",ge},{"reference",ref},{"failed",state!=0},{"half",half},{"coordinate",coordinate}};
    return *hist.emplace(name,std::move(h)).first->second;
  }
  void fill(const Hit& g,const Hit& p,const Hit* sector) {
    const int half=uint64_t(g.entry)>=split;
    std::string ref=p.kind==2?"si":p.kind==3?"s3_ring_all":"s3_sector";
    auto view=[&](const std::string& r) {
      book(g.key,r,g.failed,half,"current").Fill(g.e,p.t-g.t);
      book(g.key,r,g.failed,half,"native_ge").Fill(g.e,p.t-double(g.native));
      ++counts[std::to_string(g.key)+":"+r+":"+std::to_string(g.failed)+":"+std::to_string(half)];
    };
    view(ref);if(p.kind==3 && sector) view("s3_matched");
    if(!tree) return;
    row.ge_key=g.key;row.particle_key=p.key;row.particle_kind=p.kind;row.particle_ring=p.ring;row.particle_sector=p.sector;row.half=half;
    row.ge_entry=g.entry;row.particle_entry=p.entry;row.sector_entry=sector?sector->entry:-1;row.ge_ts=g.ts;row.particle_ts=p.ts;
    row.ge_native_ns=g.native;row.particle_native_ns=p.native;row.ge_energy_keV=g.e;row.particle_energy_keV=p.e;
    row.dt_current_ns=p.t-g.t;row.dt_native_ge_ns=p.t-double(g.native);row.dt_native_both_ns=double(p.native-g.native);
    row.particle_corrected_ns=p.t;row.ge_corrected_ns=g.t;row.ge_adc=g.adc;row.particle_adc=p.adc;
    row.ge_cfd=g.cfd;row.ge_cfds=g.cfds;row.particle_cfd=p.cfd;row.particle_cfds=p.cfds;row.ge_sr=g.sr;row.particle_sr=p.sr;
    row.ge_failed=g.failed;row.particle_failed=p.failed;row.s3_unique_match=sector!=nullptr;
    require(tree->Fill()>0,"Cannot fill pair tree");
  }
  void write() {
    file.cd();if(tree) require(tree->Write()>0,"Cannot write pairs");
    // 未测得分支也保留空谱；空谱不等于零修正。
    for(int g:GE) for(const std::string ref:{"si","s3_ring_all","s3_matched","s3_sector"})
      for(int s=0;s<2;++s) for(int h=0;h<2;++h) for(const std::string c:{"current","native_ge"}) book(g,ref,s,h,c);
    for(const auto& h:hist) require(h.second->Write()>0,"Cannot write histogram "+h.first);
  }
};
struct Engine {
  const WalkOptions& o;const Parameters& parameters;Writer& writer;
  std::vector<Hit> buffer;size_t peak_buffer=0;uint64_t ge_processed=0,matched_rings=0;
  Engine(const WalkOptions& x,const Parameters& p,Writer& w):o(x),parameters(p),writer(w) {}
  void process(double last,bool final) {
    peak_buffer=std::max(peak_buffer,buffer.size());
    require(buffer.size()<=o.max_buffer_hits,"Rolling buffer budget exceeded");
    // buffer按原生coarse时间排列；匹配索引单独按修正时间排列。
    // 唯一性要走ring->sector->另一ring两跳；相对原ring最多2个匹配窗。
    const double context=o.native_window_ns+2*SHIFT_BOUND+2*o.s3_match_ns;
    const double safe=final?std::numeric_limits<double>::infinity():last-context;
    std::vector<const Hit*> rings,sectors;
    for(const auto& h:buffer) {
      if(h.kind==3 && (parameters.group=="low" || h.ring<=8)) rings.push_back(&h);
      if(h.kind==4) sectors.push_back(&h);
    }
    auto order=[](const Hit* a,const Hit* b){return a->t<b->t || (a->t==b->t && a->entry<b->entry);};
    std::sort(rings.begin(),rings.end(),order);std::sort(sectors.begin(),sectors.end(),order);
    std::map<int64_t,const Hit*> cached_match;
    auto partners=[&](const Hit& a,const std::vector<const Hit*>& list) {
      std::vector<const Hit*> result;
      auto first=std::upper_bound(list.begin(),list.end(),a.t-o.s3_match_ns,[](double t,const Hit* h){return t<h->t;});
      for(auto it=first;it!=list.end() && (*it)->t<a.t+o.s3_match_ns;++it)
        if(std::abs((*it)->e-a.e)<o.s3_match_keV) result.push_back(*it);
      return result;
    };
    auto match=[&](const Hit& ring)->const Hit* {
      if(parameters.group!="low" && ring.ring>8) return nullptr;
      auto old=cached_match.find(ring.entry);if(old!=cached_match.end()) return old->second;
      const Hit* result=nullptr;auto s=partners(ring,sectors);
      if(s.size()==1) {auto r=partners(*s[0],rings);if(r.size()==1 && r[0]->entry==ring.entry) result=s[0];}
      cached_match[ring.entry]=result;return result;
    };
    for(auto& g:buffer) {
      if(g.kind!=1 || g.done || double(g.coarse)>safe) continue;
      const auto lo=std::upper_bound(buffer.begin(),buffer.end(),double(g.coarse)-o.native_window_ns,[](double t,const Hit& h){return t<double(h.coarse);});
      for(auto it=lo;it!=buffer.end() && double(it->coarse)<double(g.coarse)+o.native_window_ns;++it) {
        if(it->kind==1) continue;
        const Hit* sector=it->kind==3?match(*it):nullptr;
        writer.fill(g,*it,sector);if(sector) ++matched_rings;
      }
      g.done=true;++ge_processed;
    }
    if(final) {buffer.clear();return;}
    // 保留未来Ge和S3双面匹配所需的左右上下文，不在半段边界丢配对。
    const double discard=safe-context;
    auto end=std::lower_bound(buffer.begin(),buffer.end(),discard,[](const Hit& h,double t){return double(h.coarse)<t;});
    for(auto it=buffer.begin();it!=end;++it) require(it->kind!=1 || it->done,"Unprocessed Ge at buffer boundary");
    buffer.erase(buffer.begin(),end);
  }
};
void write_json_atomic(const fs::path& path,const json& value) {
  fs::path temp=path.string()+".tmp";
  {std::ofstream out(temp);require(bool(out),"Cannot create status");out<<value.dump(2)<<'\n';out.close();require(bool(out),"Cannot write status");}
  fs::rename(temp,path);
}
}

void build_ge_time_walk(const WalkOptions& input_options) {
  auto o=input_options;o.project=fs::canonical(o.project);
  require(o.ge_min>=40 && o.ge_max<=2000 && o.ge_min<o.ge_max,"Ge domain must lie in[40,2000)");
  require(o.si_min>=0 && o.si_min<o.si_max && o.s3_min>=0 && o.s3_min<o.s3_max,"Invalid particle gates");
  require(o.native_window_ns>0 && o.native_window_ns<=10000 && o.time_range_ns>0 && o.time_range_ns<=10000,"Invalid windows");
  require(o.time_bin_ns>0 && o.energy_bin_keV>0 && o.s3_match_ns>0 && o.s3_match_keV>0,"Invalid bins/matching");
  for(double n:{2*o.time_range_ns/o.time_bin_ns,(o.ge_max-o.ge_min)/o.energy_bin_keV})
    require(n>=1 && n<=10000 && std::abs(n-std::round(n))<1e-8,"Require integral number of bins<=10000");
  require(o.chunk>0 && o.chunk<=o.max_buffer_hits,"Invalid chunk/buffer limit");
  Parameters p(o.project,o.run,true);
  // ge_time_walk维护数值源；显式active保留旧流程复现，禁止隐式混用。
  fs::path profile_dir;std::string selection_source,expected_sha;
  if(!o.profile_dir.empty()) {
    require(o.profile.empty(),"Choose --profile or --profile-dir, not both");
    profile_dir=fs::canonical(o.profile_dir);o.profile=profile_dir.filename().string();
    selection_source="explicit_directory";
  } else if(o.profile.empty() || o.profile=="local") {
    const auto selected=o.project/"ge_time_walk/parameters/selected_profile.json";
    p.freeze(selected);const auto choice=json::parse(read_text(selected));
    o.profile=choice.at("profile").get<std::string>();
    require(!o.profile.empty() && o.profile.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_")==std::string::npos,"Invalid local profile");
    profile_dir=fs::canonical(selected.parent_path()/o.profile);
    expected_sha=choice.at("calibration_sha256").get<std::string>();
    selection_source="ge_time_walk/parameters/selected_profile.json";
  } else {
    if(o.profile=="active") {
      p.freeze(p.pars/"time_alignment/active.txt");
      std::istringstream active(read_text(p.pars/"time_alignment/active.txt"));std::string extra;
      require(bool(active>>o.profile) && !(active>>extra),"Invalid active profile");
      selection_source="pars/time_alignment/active.txt";
    } else selection_source="explicit_installed_version";
    require(!o.profile.empty() && o.profile.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_")==std::string::npos,"Invalid profile");
    profile_dir=fs::canonical(p.pars/"time_alignment/versions"/o.profile);
  }
  for(const auto& name:{"calibration.json","metadata.json","storage_format.txt","correction_mode.txt","source_sha256.txt"}) p.freeze(profile_dir/name);
  const auto actual_sha=p.snapshots.at((profile_dir/"calibration.json").string()).at("sha256").get<std::string>();
  std::istringstream source_hash(read_text(profile_dir/"source_sha256.txt"));std::string declared_sha,extra;
  require(bool(source_hash>>declared_sha) && !(source_hash>>extra) && declared_sha==actual_sha,"Calibration SHA does not match source_sha256.txt");
  require(expected_sha.empty() || expected_sha==actual_sha,"Local selection SHA does not match calibration.json");
  sm150::TimeAlignment timing;timing.LoadProfile(profile_dir.string(),o.run);
  require(timing.active() && timing.total(),"Require complete timing profile covering this run");
  json snapshots=p.snapshots;
  for(auto it=snapshots.begin();it!=snapshots.end();++it) if(fs::path(it.key()).parent_path()==profile_dir) {
    it.value().erase("text");it.value()["content_storage"]="versioned parameter file; SHA checked before/after processing";
  }
  json meta={{"schema","sm150_ge_walk_v1"},{"run",rn(o.run)},{"group",p.group},{"timing_profile",o.profile},
    {"timing_profile_directory",profile_dir.string()},{"timing_selection_source",selection_source},
    {"profile_sha256",p.snapshots.at((profile_dir/"calibration.json").string()).at("sha256")},
    {"parameters",p.records},{"parameter_snapshots",snapshots},{"Ge",GE},
    {"energy_gates_keV",{{"ge",{o.ge_min,o.ge_max}},{"si",{o.si_min,o.si_max}},{"s3",{o.s3_min,o.s3_max}}}},
    {"energy_gate_edges","Ge[min,max);particles(min,max)"},{"native_pair_window_ns",o.native_window_ns},
    {"pair_rule","all unique raw-entry Ge-particle pairs with strict abs(coarse_particle-coarse_Ge)<window; no nearest-time selection"},
    {"time_range_ns",o.time_range_ns},{"time_bin_ns",o.time_bin_ns},{"energy_bin_keV",o.energy_bin_keV},
    {"dt_sign","particle_minus_Ge"},{"coordinates",{{"current","particle complete time minus Ge complete time"},{"native_ge","same corrected particle minus native decoded Ge; NO Ge offset/time-walk"}}},
    {"quality","pileup=outofr=0; retain normal and failed"},{"random_jitter",false},
    {"s3_match_ns",o.s3_match_ns},{"s3_match_keV",o.s3_match_keV},
    {"s3_matching","mutually unique ring/sector within strict corrected-time and energy bounds; main matched view low ring1-24,high ring1-8; all faces also saved separately"},
    {"half_rule","Ge entry in first/second half of processed interval; sum both for final fit; boundary pairs retained"},
    {"production_active",false},{"ROOT_version",gROOT->GetVersion()}};
  if(o.inspect_parameters) {p.verify();std::cout<<meta.dump(2)<<'\n';return;}
  require(!o.input.empty() && !o.output.empty(),"Specify --input and --output");
  o.input=fs::canonical(o.input);o.output=fs::absolute(o.output).lexically_normal();
  require(!fs::exists(o.output) && !fs::is_symlink(o.output),"Output already exists");
  fs::create_directories(o.output.parent_path());
  fs::path partial=o.output.string()+".partial",progress=o.output.string()+".status.json",receipt=o.output.string()+".json";
  require(!fs::exists(partial) && !fs::is_symlink(partial) && !fs::exists(receipt),"Existing partial/receipt: inspect before retry");
  Stat before(o.input);
  std::unique_ptr<TFile> input(TFile::Open(o.input.c_str(),"READ"));require(input && !input->IsZombie(),"Cannot open raw ROOT");
  auto* tree=dynamic_cast<TTree*>(input->Get("tree"));require(tree,"Missing raw tree");
  const uint64_t total=tree->GetEntries();require(o.begin<=total,"Begin outside input");
  const uint64_t count=o.entries?std::min(o.entries,total-o.begin):total-o.begin,end=o.begin+count;
  raw reader(tree);tree->SetCacheSize(32*1024*1024);tree->AddBranchToCache("*",true);tree->StopCacheLearningPhase();
  std::unique_ptr<TFile> output(TFile::Open(partial.c_str(),"CREATE"));require(output && !output->IsZombie(),"Cannot create partial output");
  output->SetCompressionSettings(101);
  Writer writer(o,*output,o.begin+count/2);Engine engine(o,p,writer);
  Digest digest,particle_digest;std::array<InputRecord,4096> hash_records;std::array<ParticleIdentityRecord,4096> particle_records;size_t hash_used=0;
  auto flush_hash=[&](){if(hash_used){digest.add(hash_records.data(),hash_used*sizeof(InputRecord));particle_digest.add(particle_records.data(),hash_used*sizeof(ParticleIdentityRecord));hash_used=0;}};
  uint64_t rejected_quality=0,rejected_energy=0,unmapped=0,nonselected_ge=0;
  std::map<int,uint64_t> accepted;int64_t previous=std::numeric_limits<int64_t>::min();
  const auto start=std::chrono::steady_clock::now();auto last_report=start;
  auto seconds=[&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
  for(uint64_t entry=o.begin;entry<end;++entry) {
    reader.GetEntry(entry);
    hash_records[hash_used]={reader.sr,reader.pileup,reader.outofr,reader.cid,reader.sid,reader.ch,reader.evte,reader.ts,reader.cfd,reader.cfdft,reader.cfds};
    particle_records[hash_used]={reader.ts,reader.evte,reader.sr,reader.cid,reader.sid,reader.ch,uint8_t(reader.pileup),uint8_t(reader.outofr)};
    if(++hash_used==hash_records.size())flush_hash();
    require(reader.sr==100 || reader.sr==250,"Invalid sampling rate");
    require(reader.ts>=0 && reader.ts<std::numeric_limits<int64_t>::max()/10,"Invalid timestamp");
    const int64_t coarse=reader.ts*(reader.sr==250?8:10);
    require(coarse>=previous,"Raw coarse timestamps are not monotonic");previous=coarse;
    int key=10000*int(reader.cid)+100*int(reader.sid)+int(reader.ch);
    auto found=p.channels.find(key);
    if(found==p.channels.end()) ++unmapped;
    else if(reader.pileup || reader.outofr) ++rejected_quality;
    else if(found->second.kind==1 && !GE.count(key)) ++nonselected_ge;
    else {
      const auto& c=found->second;const double energy=c.energy(double(reader.evte),true);
      require(std::isfinite(energy),"Nonfinite calibrated energy");
      const bool keep=c.kind==1?energy>=o.ge_min && energy<o.ge_max:c.kind==2?energy>o.si_min && energy<o.si_max:energy>o.s3_min && energy<o.s3_max;
      if(!keep) ++rejected_energy;
      else {
        double decoded=double(coarse);
        if(c.kind==1) {require(reader.sr==250,"Ge sampling rate must be250");if(!reader.cfdft) decoded=(2.*reader.ts-reader.cfds+reader.cfd/16384.)*4.;}
        require(std::abs(decoded)<9007199254740992.,"Native timestamp exceeds exact double integer precision");
        const int64_t native=int64_t(decoded);
        const double correction=c.kind==1?timing.TotalGeOffset(key,reader.cfdft,energy):timing.TotalParticleOffsetFromAdc(key,reader.evte);
        const double corrected=double(native)+correction;
        require(std::isfinite(corrected) && std::abs(corrected-double(coarse))<=SHIFT_BOUND,"Time shift exceeds declared context bound");
        Hit h;h.key=key;h.kind=c.kind;h.ring=c.ring;h.sector=c.sector;h.entry=entry;h.coarse=coarse;h.native=native;h.ts=reader.ts;
        h.e=energy;h.t=corrected;h.adc=reader.evte;h.sr=reader.sr;h.cfd=reader.cfd;h.cfds=reader.cfds;h.failed=reader.cfdft;
        engine.buffer.push_back(h);++accepted[key];require(engine.buffer.size()<=o.max_buffer_hits,"Buffer budget exceeded");
      }
    }
    if((entry-o.begin+1)%o.chunk==0) {
      engine.process(double(coarse),false);
      if(std::chrono::duration<double>(std::chrono::steady_clock::now()-last_report).count()>10) {
        write_json_atomic(progress,{{"status","reading"},{"run",rn(o.run)},{"read_entries",entry-o.begin+1},{"requested_entries",count},{"seconds",seconds()},{"Ge_processed",engine.ge_processed}});
        std::cout<<rn(o.run)<<" "<<entry-o.begin+1<<"/"<<count<<" entries, "<<seconds()<<" s\n"<<std::flush;last_report=std::chrono::steady_clock::now();
      }
    }
  }
  engine.process(double(previous),true);require(engine.buffer.empty(),"Buffer not exhausted");
  flush_hash();
  p.verify();require(before.same(Stat(o.input)),"Input changed while reading");
  writer.write();
  meta["status"]="complete";meta["program_sha256"]=hash(read_text(fs::canonical("/proc/self/exe")));
  meta["source"]={{"path",o.input.string()},{"bytes",before.value.st_size},{"inode",before.value.st_ino},{"mtime_sec",before.value.st_mtim.tv_sec},{"mtime_nsec",before.value.st_mtim.tv_nsec},{"total_entries",total},{"entry_begin",o.begin},{"entry_end",end},{"full_run",o.begin==0 && end==total},{"read_fields_sha256",digest.finish()},{"hash_definition","sr,pileup,outofr,cid,sid,ch,evte,ts,cfd,cfdft,cfds in entry order; native little endian; no padding"}};
  meta["source"]["particle_identity_sha256"]=particle_digest.finish();
  meta["source"]["particle_identity_definition"]="packed little-endian <qHhhhhBB: ts,evte,sr,cid,sid,ch,pileup,outofr; every read record in order; matches0197 full-reader identity";
  meta["counters"]={{"read_entries",count},{"accepted",accepted},{"rejected_quality",rejected_quality},{"rejected_energy",rejected_energy},{"unmapped",unmapped},{"nonselected_Ge",nonselected_ge},{"Ge_processed",engine.ge_processed},{"matched_ring_pairs",engine.matched_rings},{"peak_buffer_hits",engine.peak_buffer}};
  meta["histograms"]=writer.descriptions;meta["pair_counts"]=writer.counts;meta["seconds"]=seconds();meta["write_pairs"]=o.write_pairs;meta["chunk_entries"]=o.chunk;
  output->cd();TObjString info(meta.dump().c_str());require(info.Write("metadata_json")>0,"Cannot write metadata");
  const int64_t saved_pair_entries=writer.tree?writer.tree->GetEntries():0;
  output->Flush();require(!output->TestBit(TFile::kWriteError),"ROOT output error");output->Close();
  {std::unique_ptr<TFile> check(TFile::Open(partial.c_str(),"READ"));require(check && !check->IsZombie(),"Output readback failed");
   for(const auto& item:writer.hist) {auto* h=dynamic_cast<TH2D*>(check->Get(item.first.c_str()));require(h && h->GetEntries()==item.second->GetEntries(),"Histogram readback mismatch");}
   if(o.write_pairs) {auto* t=dynamic_cast<TTree*>(check->Get("pairs"));require(t && t->GetEntries()==saved_pair_entries,"Pair readback mismatch");}
  }
  fs::create_hard_link(partial,o.output);fs::remove(partial);
  write_json_atomic(receipt,meta);write_json_atomic(progress,{{"status","complete"},{"run",rn(o.run)},{"read_entries",count},{"seconds",seconds()}});
  std::cout<<"COMPLETE "<<rn(o.run)<<" entries="<<count<<" Ge="<<engine.ge_processed<<" seconds="<<seconds()<<'\n';
}
