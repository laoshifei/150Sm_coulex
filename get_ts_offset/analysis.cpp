// get_ts_offset：原生 offset 测量与现用参数残差检查。
// 顺序：参数/映射 -> 原始 ADC 刻度 -> 时间 -> 排序 -> 全候选配对 -> 诊断。
// 不更新 pars，不以无数据/拟合失败的 0 代替有效 offset。
#include "analysis.h"
#include "raw.h"
#include "sm150_time_alignment.h"
#include "sm150_event_calibration.h"
#include <TFile.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TTree.h>
#include <openssl/evp.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unistd.h>
#include <sys/stat.h>
#include <vector>

namespace fs=std::filesystem;
using json=nlohmann::json;
namespace {
using namespace sm150::diagnostic;
struct Hit { double t,native,e; Long64_t entry; int key; bool failed; };
bool earlier(const Hit& a,const Hit& b) { return a.t<b.t || (a.t==b.t && a.entry<b.entry); }
struct Histograms {
  std::unique_ptr<TH1D> all, isolated, normal, failed;
  std::unique_ptr<TH2F> energy_all, energy_normal, energy_failed;
};
struct Writer {
  const OffsetOptions& o;
  TFile& file;
  int bins;
  std::map<std::string,Histograms> hist;
  json descriptions=json::object();
  TTree* pairs=nullptr;
  std::string category, name;
  int reference_key{},target_key{};
  Long64_t reference_entry{},target_entry{};
  double dt{},native_dt{},reference_energy{},target_energy{};
  ULong64_t n_reference{},n_target{};
  Bool_t isolated{},reference_failed{},target_failed{};
  Writer(const OffsetOptions& opt,TFile& f):o(opt),file(f),bins(int(std::llround(2*o.window_ns/o.bin_ns))) {
    file.cd();
    if(o.write_pairs) {
      pairs=new TTree("pairs","Each diagnostic pair; shared hits are correlated, not independent events");
      pairs->SetAutoFlush(-8*1024*1024);
      pairs->Branch("category",&category); pairs->Branch("name",&name);
      pairs->Branch("reference_key",&reference_key); pairs->Branch("target_key",&target_key);
      pairs->Branch("reference_entry",&reference_entry); pairs->Branch("target_entry",&target_entry);
      pairs->Branch("dt_ns",&dt); pairs->Branch("native_dt_ns",&native_dt);
      pairs->Branch("reference_energy_keV",&reference_energy); pairs->Branch("target_energy_keV",&target_energy);
      pairs->Branch("n_reference_candidates",&n_reference); pairs->Branch("n_target_candidates",&n_target);
      pairs->Branch("isolated",&isolated); pairs->Branch("reference_failed",&reference_failed);
      pairs->Branch("target_failed",&target_failed);
    }
  }
  Histograms& book(const std::string& id,const std::string& cat,int ref,int target,bool ge_energy=false) {
    if(hist.count(id)) return hist.at(id);
    Histograms h;
    auto h1=[&](const std::string& suffix) {
      auto p=std::make_unique<TH1D>((id+suffix).c_str(),(id+";t_{reference}-t_{target} (ns);Candidate pairs").c_str(),
                                   bins,-o.window_ns,o.window_ns);
      p->SetDirectory(nullptr); return p;
    };
    h.all=h1(""); h.isolated=h1("__isolated");
    if(ge_energy) {
      h.normal=h1("__normal"); h.failed=h1("__failed");
      auto h2=[&](const std::string& suffix) {
        auto p=std::make_unique<TH2F>((id+suffix).c_str(),(id+";Ge energy (keV);t_{Ge}-t_{particle} (ns)").c_str(),
              98,o.ge_min,o.ge_max,bins,-o.window_ns,o.window_ns);
        p->SetDirectory(nullptr); return p;
      };
      h.energy_all=h2("__energy_all"); h.energy_normal=h2("__energy_normal"); h.energy_failed=h2("__energy_failed");
    }
    descriptions[id]={{"category",cat},{"reference_key",ref},{"target_key",target},
      {"dt_definition","reference_minus_target"},{"isolated_definition","one candidate in each direction for this channel pair"},
      {"ge_energy_diagnostics",ge_energy}};
    return hist.emplace(id,std::move(h)).first->second;
  }
  void fill(const std::string& id,const std::string& cat,const Hit& a,const Hit& b,
            size_t nr,size_t nt,bool ge_energy=false) {
    auto& h=book(id,cat,a.key,ge_energy?-1:b.key,ge_energy);
    const double d=a.t-b.t; const bool iso=nr==1 && nt==1;
    h.all->Fill(d); if(iso) h.isolated->Fill(d);
    if(ge_energy) {
      (a.failed?h.failed:h.normal)->Fill(d); h.energy_all->Fill(a.e,d);
      (a.failed?h.energy_failed:h.energy_normal)->Fill(a.e,d);
    }
    if(pairs) {
      category=cat;name=id;reference_key=a.key;target_key=b.key;
      reference_entry=a.entry;target_entry=b.entry;dt=d;native_dt=a.native-b.native;
      reference_energy=a.e;target_energy=b.e;n_reference=nr;n_target=nt;
      isolated=iso;reference_failed=a.failed;target_failed=b.failed;
      require(pairs->Fill()>0,"Pair tree write failed");
    }
  }
  void write() {
    file.cd(); if(pairs) require(pairs->Write()>0,"Cannot write pair tree");
    for(const auto& kv:hist) {
      const auto& h=kv.second;
      for(auto* p:{h.all.get(),h.isolated.get(),h.normal.get(),h.failed.get()})
        if(p) require(p->Write()>0,"Cannot write histogram");
      for(auto* p:{h.energy_all.get(),h.energy_normal.get(),h.energy_failed.get()})
        if(p) require(p->Write()>0,"Cannot write energy histogram");
    }
  }
};
std::string legacy_name(int key) {
  std::ostringstream s; s<<(key<10000?"ge":"si")<<"_sid"<<(key%10000/100)<<"_ch"<<std::setw(2)<<std::setfill('0')<<key%100;
  return s.str();
}
// 分块只影响内存，不改变配对。保留过去 window 和未来 2*window，才能判断双向多重性。
// 1e6ns 是经逐事件断言的最大时间移动界，不是符合窗；与 build 使用同一界限。
struct PairEngine {
  const OffsetOptions& o; const Parameters& p; Writer& w;
  std::vector<Hit> buffer;
  double consumed=-std::numeric_limits<double>::infinity();
  size_t peak_buffer=0;
  static constexpr double shift_bound=1000000.;
  PairEngine(const OffsetOptions& x,const Parameters& y,Writer& z):o(x),p(y),w(z) {}
  void process(double last_coarse,bool final) {
    peak_buffer=std::max(peak_buffer,buffer.size());
    require(buffer.size()<=o.max_buffer_hits,"Rolling buffer budget exceeded; reduce --chunk or inspect time ordering");
    std::sort(buffer.begin(),buffer.end(),earlier);
    std::map<int,std::vector<double>> times,normal;
    for(const auto& h:buffer) { times[h.key].push_back(h.t); if(!h.failed) normal[h.key].push_back(h.t); }
    auto count=[&](int key,double center,bool only_normal=false) {
      const auto& v=only_normal?normal[key]:times[key];
      auto lo=std::upper_bound(v.begin(),v.end(),center-o.window_ns);
      auto hi=std::lower_bound(v.begin(),v.end(),center+o.window_ns);
      return size_t(hi-lo);
    };
    const double safe=final?std::numeric_limits<double>::infinity():last_coarse-shift_bound-2*o.window_ns;
    require(final || safe>=consumed,"Nonmonotonic processing watermark");
    auto emit=[&](const std::string& id,const std::string& cat,const Hit& ref,const Hit& target,
                  bool ref_normal=false,bool ge_energy=false) {
      if(ref_normal && ref.failed) return;
      // n_reference = number of allowed reference hits near this target.
      w.fill(id,cat,ref,target,count(ref.key,target.t,ref_normal),count(target.key,ref.t),ge_energy);
    };
    auto si_ref=[](int k){return k==10400 || k==10600;};
    auto s3_ref=[&](int k) {
      return k==11102 || k==11200 || k==(p.group=="low"?10200:10900) || k==(p.group=="low"?10300:11000);
    };
    for(size_t i=0;i<buffer.size();++i) {
      const auto& a=buffer[i];
      if(a.t<=consumed) continue;
      if(a.t>safe) break;
      for(size_t j=i+1;j<buffer.size() && buffer[j].t-a.t<o.window_ns;++j) {
        const auto& b=buffer[j]; if(a.key==b.key) continue;
        const int ka=p.channels.at(a.key).kind,kb=p.channels.at(b.key).kind;
        if(o.particle_network) {
          // Si作为参考的独立诊断：不要求Ge，不让同一对的两个方向重复填谱。
          if(ka==2 && kb==2) {
            const auto& r=a.key<b.key?a:b; const auto& s=a.key<b.key?b:a;
            emit("si_si_"+std::to_string(r.key)+"_"+std::to_string(s.key),"si_si",r,s);
          } else if((ka==2 && (kb==3 || kb==4)) || (kb==2 && (ka==3 || ka==4))) {
            const auto& r=ka==2?a:b; const auto& s=ka==2?b:a;
            emit("si_s3_"+std::to_string(r.key)+"_"+std::to_string(s.key),"si_s3",r,s);
          }
          continue;
        }
        if(o.particle_calibration) {
          // 单独的标定诊断模式：保留所有粒子通道，用既有20个Ge作参考。
          // normal/failed及合计同时保存；不改变普通诊断模式和生产build。
          static const std::set<int> refs={204,205,206,209,211,212,213,214,215,
            304,305,306,308,310,311,506,507,508,509,510};
          if((ka==1)!=(kb==1)) {
            const auto& g=ka==1?a:b; const auto& q=ka==1?b:a;
            if(refs.count(g.key)) {
              std::string id="particle_"+std::to_string(q.key)+"_ge";
              emit(id,"particle_calibration_ge",g,q);
              emit(id+(g.failed?"_failed":"_normal"),"particle_calibration_ge",g,q);
            }
          } else if(((ka==3 && kb==4)||(ka==4 && kb==3)) && std::abs(a.e-b.e)<2500.) {
            const auto& r=ka==3?a:b; const auto& s=ka==4?a:b;
            emit("s3_rs_"+std::to_string(r.key)+"_"+std::to_string(s.key),"s3_energy_matched",r,s);
          }
          continue;
        }
        if(ka==1 || kb==1) {
          if(ka==1 && kb==1) {
            if(a.key==o.ge_reference) emit(legacy_name(b.key),"ge_ge",a,b,true);
            else if(b.key==o.ge_reference) emit(legacy_name(a.key),"ge_ge",b,a,true);
          } else {
            const auto& g=ka==1?a:b; const auto& q=ka==1?b:a;
            int kind=p.channels.at(q.key).kind;
            std::string family=kind==2?"si":kind==3?"s3_ring":"s3_sector";
            emit("ge"+std::to_string(g.key)+"_"+family,"ge_particle",g,q,false,true);
            if(g.key==o.ge_reference) emit(legacy_name(q.key),"particle_ge",g,q,true);
          }
        } else if(ka==2 && kb==2) {
          if(si_ref(a.key) || si_ref(b.key)) {
            const bool ar=si_ref(a.key) && (!si_ref(b.key) || a.key<b.key);
            const auto& r=ar?a:b; const auto& t=ar?b:a;
            emit("si_si_"+std::to_string(r.key)+"_"+std::to_string(t.key),"si_si",r,t);
          }
        } else if((ka==3 && kb==4)||(ka==4 && kb==3)) {
          const auto& r=ka==3?a:b; const auto& s=ka==4?a:b;
          emit("s3_rs_"+std::to_string(r.key)+"_"+std::to_string(s.key),"s3_ring_sector",r,s);
        } else if(ka==kb && (ka==3 || ka==4) && (s3_ref(a.key)||s3_ref(b.key))) {
          const bool ar=s3_ref(a.key) && (!s3_ref(b.key) || a.key<b.key);
          const auto& r=ar?a:b; const auto& t=ar?b:a;
          emit("s3_same_"+std::to_string(r.key)+"_"+std::to_string(t.key),"s3_same_face",r,t);
        } else if(ka==2 || kb==2) {
          const auto& s=ka==2?a:b; const auto& q=ka==2?b:a;
          if(si_ref(s.key) || s3_ref(q.key))
            emit("si_s3_"+std::to_string(s.key)+"_"+std::to_string(q.key),"si_s3",s,q);
        }
      }
    }
    consumed=safe;
    // 未来未处理 hit 的反向候选仍可能在 safe-window 之后，必须保留。
    auto end=std::upper_bound(buffer.begin(),buffer.end(),safe-o.window_ns,
                             [](double t,const Hit& h){return t<h.t;});
    buffer.erase(buffer.begin(),end);
  }
};
// 全统计粒子模式：单遍遍历所有原始条目，在线累积整数直方图。
// 粒子原生时间只有coarse项，单调排序允许用两个滚动队列枚举全部候选。
// 每对只在较后原始条目到达时计数；跨前/后半边界的配对归入较后半。
struct FullCounts {
  size_t rows;
  std::vector<uint64_t> bins;
  FullCounts(size_t n):rows(n),bins(2*n*500,0) {}
  void fill(int half,size_t row,int bin) { ++bins[(size_t(half)*rows+row)*500+size_t(bin)]; }
  void write(const std::string& name,TFile& file) const {
    file.cd();
    for(int half=0;half<2;++half) {
      const std::string id=name+"_"+std::to_string(half);
      TH2D h(id.c_str(),(name+";time difference (ns);row index").c_str(),500,-1000,1000,int(rows),0,rows);
      h.SetDirectory(nullptr); uint64_t sum=0;
      for(size_t row=0;row<rows;++row) for(int bin=0;bin<500;++bin) {
        uint64_t v=bins[(size_t(half)*rows+row)*500+size_t(bin)];
        require(v<9007199254740992ULL,"Count exceeds exact ROOT double range");
        h.SetBinContent(bin+1,int(row)+1,double(v));sum+=v;
      }
      h.SetEntries(double(sum));require(h.Write()>0,"Cannot write full-statistics histogram");
    }
  }
  std::array<uint64_t,2> totals() const {
    std::array<uint64_t,2> out{};
    for(size_t i=0;i<bins.size();++i) out[i/(rows*500)]+=bins[i];
    return out;
  }
};
struct FullHit { int64_t t; int index; double energy,offset; int ebin; };
void full_particle(const OffsetOptions& o,Parameters& p) {
  require(o.mode=="raw" && o.profile.empty() && o.energy_stage=="final",
    "particle-full fixes native selection and the documented Si/S3 baseline; use raw/final without profile");
  require(!o.particle_network && !o.particle_calibration && !o.write_pairs,"particle-full has its own compact histogram output");
  require(o.window_ns==2000 && o.bin_ns==4 && o.si_min==2500 && o.s3_min==2500 &&
          o.si_max==100000 && o.s3_max==100000,"particle-full requires window2000/bin4/particle2.5-100MeV");
  require(!o.input.empty() && !o.output.empty(),"Specify input and output");
  const fs::path input_path=fs::canonical(o.input), output_path=fs::absolute(o.output).lexically_normal();
  const fs::path partial=output_path.string()+".partial";
  require(!fs::exists(output_path) && !fs::is_symlink(output_path) && !fs::exists(partial) && !fs::is_symlink(partial),"Output/partial exists");
  fs::create_directories(output_path.parent_path());
  const fs::path base=p.pars/"time_alignment/versions/high_low_beam_failed_centers_20261021/calibration.json";
  const fs::path s3dir=p.pars/"time_alignment/versions/particle_channels_20261023";
  p.freeze(base);
  for(const auto& f:fs::recursive_directory_iterator(s3dir)) if(f.is_regular_file()) p.freeze(f.path());
  const auto constants=json::parse(read_text(base)).at("channel_offsets").at(p.group);
  sm150::TimeAlignment s3;s3.LoadProfile(s3dir.string(),o.run);require(s3.total(),"Missing complete S3 baseline");
  std::vector<int> keys,si,s3keys,rings,sectors;
  std::array<int,16384> lookup;lookup.fill(-1);
  for(const auto& kv:p.channels) if(kv.second.kind>=2) {
    lookup.at(kv.first)=int(keys.size());keys.push_back(kv.first);
    if(kv.second.kind==2) si.push_back(kv.first);
    else {s3keys.push_back(kv.first);(kv.second.kind==3?rings:sectors).push_back(kv.first);}
  }
  require(si.size()==96 && s3keys.size()==56 && rings.size()==24 && sectors.size()==32,"Unexpected particle mapping");
  std::vector<int> si_index(keys.size(),-1),s3_index(keys.size(),-1),face(keys.size(),-1),face_index(keys.size(),-1);
  for(size_t i=0;i<si.size();++i) si_index[lookup[si[i]]]=int(i);
  for(size_t i=0;i<s3keys.size();++i) s3_index[lookup[s3keys[i]]]=int(i);
  for(size_t i=0;i<rings.size();++i) {face[lookup[rings[i]]]=0;face_index[lookup[rings[i]]]=int(i);}
  for(size_t i=0;i<sectors.size();++i) {face[lookup[sectors[i]]]=1;face_index[lookup[sectors[i]]]=int(i);}
  const int reference=p.group=="low"?11007:p.group=="high1"?10209:10200;
  const std::array<double,9> energy_edges={2500,5000,8000,12000,15000,25000,40000,60000,100000};
  std::vector<double> energies(keys.size()*65536),offsets(energies.size());
  for(size_t i=0;i<keys.size();++i) {
    const auto& c=p.channels.at(keys[i]);
    for(size_t adc=0;adc<65536;++adc) {
      const size_t at=i*65536+adc;const double e=c.energy(adc,true);energies[at]=e;
      if(e>2500 && e<100000) {
        offsets[at]=c.kind==2?constants.at(std::to_string(keys[i])).get<double>():s3.TotalParticleOffsetFromAdc(keys[i],adc);
        require(std::isfinite(offsets[at]) && std::abs(offsets[at])<1000000,"Invalid baseline offset");
      }
    }
  }
  FullCounts ss(4560),ss_high(4560),cross(96*56),cross_high(96*56),energy(96*2*8),cross_energy(96*2*8*8),ref_energy(96*8*8),rs(24*32),rs_high(24*32);
  std::array<std::array<int,96>,96> pair_index{};int kpair=0;
  for(int i=0;i<96;++i) for(int j=i+1;j<96;++j) pair_index[i][j]=kpair++;
  Stat before(input_path);
  std::unique_ptr<TFile> input(TFile::Open(input_path.c_str(),"READ"));require(input && !input->IsZombie(),"Cannot read source");
  auto* tree=dynamic_cast<TTree*>(input->Get("tree"));require(tree,"Missing source tree");
  const uint64_t total=tree->GetEntries();require(o.begin<total,"Begin outside source");
  const uint64_t end=o.entries?o.begin+std::min(o.entries,total-o.begin):total;
  // 对读取的八个标量字段逐条打包后批量SHA；不读取或依赖Ge CFD字段。
#pragma pack(push,1)
  struct Record { Long64_t ts;UShort_t evte;Short_t sr,cid,sid,ch;Bool_t pileup,outofr; } record{};
#pragma pack(pop)
  static_assert(sizeof(Record)==20,"Unexpected packed record layout");
  tree->SetBranchStatus("*",0);
  auto bind=[&](const char* name,const char* type,void* ptr) {
    auto* leaf=tree->GetLeaf(name);require(tree->GetBranch(name) && leaf && !leaf->GetLeafCount() &&
      leaf->GetLenStatic()==1 && std::string(leaf->GetTypeName())==type,"Wrong source scalar field");
    tree->SetBranchStatus(name,1);require(tree->SetBranchAddress(name,ptr)>=0,"Cannot bind source field");
  };
  bind("ts","Long64_t",&record.ts);bind("evte","UShort_t",&record.evte);
  bind("sr","Short_t",&record.sr);bind("cid","Short_t",&record.cid);bind("sid","Short_t",&record.sid);bind("ch","Short_t",&record.ch);
  bind("pileup","Bool_t",&record.pileup);bind("outofr","Bool_t",&record.outofr);
  // 上述八项中的原始ADC、时间及通道/质量全部保留；采用完整标量记录哈希。
  tree->SetCacheSize(64*1024*1024);
  for(const auto* name:{"ts","evte","sr","cid","sid","ch","pileup","outofr"}) tree->AddBranchToCache(name,true);
  tree->StopCacheLearningPhase();
  std::unique_ptr<TFile> output(TFile::Open(partial.c_str(),"CREATE"));require(output && !output->IsZombie(),"Cannot create full output");
  output->SetCompressionSettings(504);
  std::vector<Record> hash_buffer(65536);size_t hash_used=0;Digest digest;
  std::deque<FullHit> sibuffer,s3buffer;size_t peak_buffer=0;
  uint64_t read=0,quality=0,unmapped=0,ge_skipped=0,energy_rejected=0;
  std::array<uint64_t,2> pair_ss{},pair_cross{},pair_rs{},read_half{};
  std::vector<uint64_t> accepted(keys.size(),0),singles(keys.size()*200,0);
  int64_t previous=-1;
  const auto start=std::chrono::steady_clock::now();double last_progress=-100;
  auto progress=[&](const std::string& status) {
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    json j={{"status",status},{"run",rn(o.run)},{"pid",::getpid()},{"input",input_path.string()},
      {"entry_begin",o.begin},{"entry_end",end},{"entries_read",read},{"total_entries",total},
      {"fraction",double(read)/(end-o.begin)},{"seconds",seconds},{"records_per_second",seconds>0?read/seconds:0},
      {"si_si_pairs",pair_ss},{"si_s3_pairs",pair_cross},{"s3_rs_pairs",pair_rs}};
    fs::path path=output_path.string()+".progress.json",tmp=path.string()+".tmp";
    {std::ofstream f(tmp);require(bool(f),"Cannot write progress");f<<j.dump(2)<<'\n';require(bool(f),"Progress write failed");}
    fs::rename(tmp,path);last_progress=seconds;
  };
  auto emit_ss=[&](const FullHit& u,const FullHit& v,int half) {
    const auto& a=si_index[u.index]<si_index[v.index]?u:v;const auto& b=si_index[u.index]<si_index[v.index]?v:u;
    const int ai=si_index[a.index],bi=si_index[b.index];if(ai==bi)return;
    ++pair_ss[half];const double dt=double(a.t-b.t)+a.offset-b.offset;
    if(dt< -1000 || dt>=1000)return;
    const int bin=int(std::floor((dt+1000)/4));
    const size_t row=pair_index[ai][bi];ss.fill(half,row,bin);
    if(a.energy>15000 && b.energy>15000)ss_high.fill(half,row,bin);
    if(keys[a.index]==reference || keys[b.index]==reference) {
      const auto& ref=keys[a.index]==reference?a:b;const auto& target=keys[a.index]==reference?b:a;
      const double d=double(target.t-ref.t)+target.offset-ref.offset;
      if(d>=-1000 && d<1000)ref_energy.fill(half,(si_index[target.index]*8+target.ebin)*8+ref.ebin,int(std::floor((d+1000)/4)));
    }
  };
  auto emit_cross=[&](const FullHit& a,const FullHit& b,int half) {
    ++pair_cross[half];const double dt=double(a.t-b.t)+a.offset-b.offset;
    if(dt< -1000 || dt>=1000)return;
    const int bin=int(std::floor((dt+1000)/4));
    const size_t row=si_index[a.index]*56+s3_index[b.index];cross.fill(half,row,bin);
    if(a.energy>15000 && b.energy>25000)cross_high.fill(half,row,bin);
    const size_t erow=(si_index[a.index]*2+face[b.index])*8+a.ebin;
    energy.fill(half,erow,bin);cross_energy.fill(half,erow*8+b.ebin,bin);
  };
  auto emit_rs=[&](const FullHit& u,const FullHit& v,int half) {
    if(face[u.index]==face[v.index] || std::abs(u.energy-v.energy)>=2500)return;
    const auto& a=face[u.index]==0?u:v;const auto& b=face[u.index]==0?v:u;
    ++pair_rs[half];const double dt=double(a.t-b.t)+a.offset-b.offset;
    if(dt< -1000 || dt>=1000)return;
    const int bin=int(std::floor((dt+1000)/4));
    const size_t row=face_index[a.index]*32+face_index[b.index];rs.fill(half,row,bin);
    if(a.energy>25000 && b.energy>25000)rs_high.fill(half,row,bin);
  };
  progress("reading");
  for(uint64_t entry=o.begin;entry<end;++entry) {
    if(tree->GetEntry(entry)<=0)throw std::runtime_error("Source read failed at "+std::to_string(entry));
    hash_buffer[hash_used++]=record;if(hash_used==hash_buffer.size()){digest.add(hash_buffer.data(),hash_used*sizeof(Record));hash_used=0;}
    if((record.sr!=100 && record.sr!=250) || record.cid<0 || record.cid>1 || record.sid<2 || record.sid>12 || record.ch<0 || record.ch>=16 || record.ts<0 || record.ts>std::numeric_limits<int64_t>::max()/10)
      throw std::runtime_error("Invalid source record at "+std::to_string(entry));
    const int64_t coarse=record.ts*(record.sr==250?8:10);
    if(coarse<previous)throw std::runtime_error("Coarse time decreased at "+std::to_string(entry));
    previous=coarse;
    const int half=entry>=total/2; ++read_half[half];++read;
    const int key=record.cid*10000+record.sid*100+record.ch;
    if(record.cid==0)++ge_skipped;
    else if(lookup[key]<0)++unmapped;
    else if(record.pileup || record.outofr)++quality;
    else {
      const int i=lookup[key];const size_t at=size_t(i)*65536+record.evte;const double e=energies[at];
      if(e>=0 && e<100000)++singles[size_t(i)*200+size_t(e/500)];
      if(!(e>2500 && e<100000))++energy_rejected;
      else {
        ++accepted[i];FullHit hit{coarse,i,e,offsets[at],int(std::upper_bound(energy_edges.begin(),energy_edges.end(),e)-energy_edges.begin()-1)};
        while(!sibuffer.empty() && sibuffer.front().t<=coarse-2000)sibuffer.pop_front();
        while(!s3buffer.empty() && s3buffer.front().t<=coarse-2000)s3buffer.pop_front();
        if(si_index[i]>=0) {
          for(const auto& b:sibuffer)emit_ss(hit,b,half);
          for(const auto& b:s3buffer)emit_cross(hit,b,half);
          sibuffer.push_back(hit);
        } else {
          for(const auto& a:sibuffer)emit_cross(a,hit,half);
          for(const auto& b:s3buffer)emit_rs(hit,b,half);
          s3buffer.push_back(hit);
        }
        peak_buffer=std::max(peak_buffer,sibuffer.size()+s3buffer.size());
        if(peak_buffer>o.max_buffer_hits)throw std::runtime_error("Particle rolling buffer exceeded");
      }
    }
    if(read%1000000==0 && std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()-last_progress>10)progress("reading");
  }
  if(hash_used)digest.add(hash_buffer.data(),hash_used*sizeof(Record));
  require(before.same(Stat(input_path)),"Source changed during complete read");p.verify();progress("writing");
  json hist_counts=json::object();
  for(const auto& kv:std::vector<std::pair<std::string,FullCounts*>>{{"ss",&ss},{"ss_high",&ss_high},{"cross",&cross},{"cross_high",&cross_high},{"energy",&energy},{"cross_energy",&cross_energy},{"reference_energy",&ref_energy},{"rs",&rs},{"rs_high",&rs_high}}) {
    kv.second->write(kv.first,*output);hist_counts[kv.first]=kv.second->totals();
  }
  json metadata={{"schema","sm150_particle_full_v1"},{"status","complete"},{"run",rn(o.run)},{"group",p.group},
    {"source",{{"path",input_path.string()},{"bytes",before.value.st_size},{"mtime_sec",before.value.st_mtim.tv_sec},{"mtime_nsec",before.value.st_mtim.tv_nsec},{"total_entries",total},{"entry_begin",o.begin},{"entry_end",end},{"full_run",o.begin==0 && end==total},{"read_fields_sha256",digest.finish()},
    {"hash_definition","packed little-endian <qHhhhhBB: ts,evte,sr,cid,sid,ch,pileup,outofr; every read record in order"}}},
    {"read_entries",read},{"read_half",read_half},{"split_rule","first/second entry half; each pair assigned to later entry; no records discarded; fit final sums both"},
    {"particle_keys",keys},{"si_keys",si},{"s3_keys",s3keys},{"s3_ring_keys",rings},{"s3_sector_keys",sectors},{"reference_key",reference},
    {"histogram_counts",hist_counts},{"histogram_bin_ns",4},{"histogram_time_range_ns",{-1000,1000}},{"native_pair_window_ns",2000},
    {"energy_edges_keV",energy_edges},{"parameters",p.records},{"parameter_snapshots",p.snapshots},
    {"Si_baseline","legacy channel constants high_low_beam_failed_centers_20261021"},{"S3_baseline","particle_channels_20261023 complete internal function"},
    {"Ge_used",false},{"quality","pileup=outofr=0; particle energy strict2.5-100MeV"},
    {"counters",{{"ge_skipped",ge_skipped},{"unmapped",unmapped},{"rejected_quality",quality},{"rejected_energy",energy_rejected},{"accepted",accepted},{"peak_buffer_hits",peak_buffer},{"si_si_pairs",pair_ss},{"si_s3_pairs",pair_cross},{"s3_rs_pairs",pair_rs}}},
    {"singles_500keV_bins",singles},{"program_sha256",hash(read_text(fs::canonical("/proc/self/exe")))},{"ROOT_version",gROOT->GetVersion()},
    {"seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()},{"candidate_status","analysis_only_not_adopted"}};
  TObjString record_json(metadata.dump().c_str());require(record_json.Write("metadata_json")>0,"Cannot write metadata");
  output->Flush();require(!output->TestBit(TFile::kWriteError),"Full output write error");output->Close();
  {std::unique_ptr<TFile> check(TFile::Open(partial.c_str(),"READ"));require(check && !check->IsZombie(),"Cannot read completed output");
    auto* m=dynamic_cast<TObjString*>(check->Get("metadata_json"));require(m && json::parse(m->GetString().Data()).at("read_entries")==read,"Metadata readback failed");
    for(auto it=hist_counts.begin();it!=hist_counts.end();++it) for(int h=0;h<2;++h) {
      auto* hist=dynamic_cast<TH2D*>(check->Get((it.key()+"_"+std::to_string(h)).c_str()));
      require(hist && hist->Integral()==it.value()[h].get<double>(),"Histogram count readback failed");
    }
  }
  fs::create_hard_link(partial,output_path);fs::remove(partial);progress("complete");
  std::cout<<json({{"status","complete"},{"output",output_path.string()},{"full_run",o.begin==0 && end==total},{"read_entries",read},{"seconds",metadata["seconds"]}}).dump()<<'\n';
}
} // namespace

void analysis::Process() {
  auto o=options_;
  require(o.mode=="raw" || o.mode=="residual","Mode must be raw or residual");
  require(o.profile.empty() || o.mode=="residual","An explicit profile requires residual mode");
  require(!(o.particle_network && o.particle_calibration),"Choose one compact diagnostic mode");
  require(o.energy_stage=="base" || o.energy_stage=="final","Energy stage must be base or final");
  require(o.window_ns>0 && o.window_ns<=100000 && o.bin_ns>0 && o.bin_ns<=o.window_ns,"Invalid time window/bin");
  const double bins=2*o.window_ns/o.bin_ns;
  require(bins<=10000 && std::abs(bins-std::round(bins))<1e-8,"Window must contain an integer number of bins (<=10000)");
  require(o.chunk>0 && o.max_buffer_hits>0 && o.chunk<=o.max_buffer_hits,"Invalid chunk/buffer limit");
  require(o.ge_min>=0 && o.ge_min<o.ge_max && o.si_min>=0 && o.si_min<o.si_max &&
          o.s3_min>=0 && o.s3_min<o.s3_max,"Invalid energy gates");
  o.project=fs::canonical(o.project);
  Parameters p(o.project,o.run,o.energy_stage=="final");
  if(o.particle_full && !o.inspect_parameters) { full_particle(o,p);return; }
  require(p.channels.count(o.ge_reference) && p.channels.at(o.ge_reference).kind==1,"Invalid Ge reference");
  sm150::TimeAlignment timing;
  if(o.mode=="residual") {
    require(o.energy_stage=="final","Residual mode requires final energy, consistent with active timing profile");
    p.freeze(p.pars/"time_alignment/active.txt");
    std::istringstream active(p.snapshots.at((p.pars/"time_alignment/active.txt").string()).at("text").get<std::string>());
    std::string profile; active>>profile;
    if(!o.profile.empty()) profile=o.profile;
    require(!profile.empty() && profile.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_")==std::string::npos,"Invalid active profile");
    fs::path dir=p.pars/"time_alignment/versions"/profile;
    for(const auto& item:fs::recursive_directory_iterator(dir)) if(item.is_regular_file()) p.freeze(item.path());
    timing.LoadProfile(dir.string(),o.run);
    require(timing.active() && timing.total(),"Residual mode requires a complete current profile");
  }
  json metadata={{"schema","sm150_offset_v2"},{"run",rn(o.run)},{"group",p.group},{"mode",o.mode},
    {"energy_stage",o.energy_stage},{"time_definition","int64(native CFD/coarse ns), matching build; optional complete offset added once"},
    {"window_ns",o.window_ns},{"bin_ns",o.bin_ns},{"window_rule","abs(dt)<window"},{"ge_reference",o.ge_reference},
    {"reference_cfd","normal_only_for_legacy_reference_spectra"},{"dt_sign","reference_minus_target"},
    {"energy_gates_keV",{{"ge",{o.ge_min,o.ge_max}},{"si",{o.si_min,o.si_max}},{"s3",{o.s3_min,o.s3_max}}}},
    {"energy_gate_edges","Ge [min,max); particles (min,max)"},{"quality","reject pileup and outofr; keep normal and failed"},
    {"particle_refs",{{"si",{10400,10600}},{"s3_ring",{11102,11200}},
                     {"s3_sector",p.group=="low"?json::array({10200,10300}):json::array({10900,11000})}}},
    {"parameters",p.records},{"parameter_snapshots",p.snapshots},{"timing_profile",timing.profile()},
    {"random_jitter",false},{"ROOT_version",gROOT->GetVersion()},
    {"particle_calibration",o.particle_calibration},
    {"particle_network",o.particle_network},
    {"particle_calibration_rules","Ge reference keys:204,205,206,209,211,212,213,214,215,304,305,306,308,310,311,506,507,508,509,510; S3 |E_ring-E_sector|<2500keV; all candidates retained"},
    {"candidate_status","diagnostic_only_not_adopted"},{"notes","All pairs can share hits; isolated means mutually unique within one channel pair, not whole-event uniqueness."}};
  if(o.inspect_parameters) { p.verify(); std::cout<<metadata.dump(2)<<'\n'; return; }
  metadata["program_sha256"]=hash(read_text(fs::canonical("/proc/self/exe")));
  require(!o.input.empty() && !o.output.empty(),"Specify --input and --output");
  o.input=fs::canonical(o.input); o.output=fs::absolute(o.output).lexically_normal();
  require(!fs::exists(o.output) && !fs::is_symlink(o.output),"Output already exists");
  fs::create_directories(o.output.parent_path());
  const fs::path partial=o.output.string()+".partial";
  require(!fs::exists(partial) && !fs::is_symlink(partial),"Partial output exists; inspect it before retrying");
  Stat before(o.input);
  std::unique_ptr<TFile> input(TFile::Open(o.input.c_str(),"READ"));
  require(input && !input->IsZombie(),"Cannot open input ROOT");
  auto* tree=dynamic_cast<TTree*>(input->Get("tree")); require(tree,"Missing input tree");
  const uint64_t total=tree->GetEntries();
  require(o.begin<total,"Begin entry outside input");
  const uint64_t count=o.entries?std::min(o.entries,total-o.begin):total-o.begin;
  const uint64_t end=o.begin+count;
  raw rd(tree);
  // CREATE 的 O_EXCL 语义防止并发同名输出；正式文件只在全部检查后发布。
  std::unique_ptr<TFile> output(TFile::Open(partial.c_str(),"CREATE"));
  require(output && !output->IsZombie(),"Cannot create partial output");
  Writer writer(o,*output);
  PairEngine engine(o,p,writer);
  Digest digest;
  std::map<int,uint64_t> accepted,normal_counts,failed_counts;
  uint64_t rejected_quality=0,rejected_energy=0,unmapped=0;
  int64_t previous=std::numeric_limits<int64_t>::min();
  const auto start=std::chrono::steady_clock::now();
  for(uint64_t entry=o.begin;entry<end;++entry) {
    rd.GetEntry(entry);
    // 对实际读入的原始标量字段求 SHA256；不是未读取整份文件的 SHA256。
    digest.add(rd.sr);digest.add(rd.pileup);digest.add(rd.outofr);digest.add(rd.cid);digest.add(rd.sid);
    digest.add(rd.ch);digest.add(rd.evte);digest.add(rd.ts);digest.add(rd.cfd);digest.add(rd.cfdft);digest.add(rd.cfds);
    require(rd.sr==100 || rd.sr==250,"Unsupported sampling rate");
    require((rd.cid==0 || rd.cid==1) && rd.sid>=2 && rd.sid<=12 && rd.ch>=0 && rd.ch<16,
            "Invalid hardware channel at entry "+std::to_string(entry));
    require(rd.ts>=0 && rd.ts<=std::numeric_limits<int64_t>::max()/10,"Invalid coarse timestamp");
    const int64_t coarse=rd.ts*(rd.sr==250?8:10);
    require(coarse>=previous,"Coarse timestamp order decreased at entry "+std::to_string(entry));
    previous=coarse;
    const int key=10000*rd.cid+100*rd.sid+rd.ch;
    auto found=p.channels.find(key);
    if(found==p.channels.end()) ++unmapped;
    else if(rd.pileup || rd.outofr) ++rejected_quality;
    else {
      const auto& c=found->second;
      const double energy=c.energy(rd.evte,o.energy_stage=="final");
      require(std::isfinite(energy),"Nonfinite energy");
      bool keep=c.kind==1 ? energy>=o.ge_min && energy<o.ge_max :
                c.kind==2 ? energy>o.si_min && energy<o.si_max : energy>o.s3_min && energy<o.s3_max;
      if(o.particle_network && c.kind==1) keep=false;
      if(!keep) ++rejected_energy;
      else {
        double native=double(coarse);
        if(c.kind==1) {
          require(rd.sr==250,"Ge sample rate must be 250");
          native=rd.cfdft?double(coarse):(2.*rd.ts-rd.cfds+rd.cfd/16384.)*4.;
        }
        require(std::abs(native)<9007199254740992.,"Timestamp exceeds exact double integer precision");
        native=double(int64_t(native)); // 保持当前 build 的原生时间量化定义。
        const double correction=o.mode=="raw"?0.:c.kind==1 ?
            timing.TotalGeOffset(key,rd.cfdft,energy):timing.TotalParticleOffsetFromAdc(key,rd.evte);
        const double t=native+correction;
        require(std::isfinite(t) && std::abs(t-coarse)<=PairEngine::shift_bound,"Timing shift exceeds streaming bound");
        engine.buffer.push_back({t,native,energy,Long64_t(entry),key,bool(rd.cfdft)});
        require(engine.buffer.size()<=o.max_buffer_hits,"Rolling buffer budget exceeded");
        ++accepted[key]; ++(rd.cfdft?failed_counts:normal_counts)[key];
      }
    }
    if((entry-o.begin+1)%o.chunk==0) engine.process(double(coarse),false);
  }
  engine.process(double(previous),true);
  require(engine.buffer.empty(),"Final buffer not consumed");
  // 保留所有物理通道的空谱及显式状态，参考缺失不阻止独立粒子诊断。
  json coverage=json::object();
  for(const auto& kv:p.channels) {
    int k=kv.first;
    writer.book(legacy_name(k),kv.second.kind==1?"ge_ge":"particle_ge",o.ge_reference,k);
    coverage[std::to_string(k)]={{"accepted",accepted[k]},{"normal",normal_counts[k]},{"failed",failed_counts[k]},
      {"status",k==o.ge_reference?(normal_counts[k]?"reference":"reference_missing"):
        !normal_counts[o.ge_reference]?"reference_missing":!accepted[k]?"no_selected_hits":"measured_histogram"}};
  }
  require(before.same(Stat(o.input)),"Input changed while reading"); p.verify();
  metadata["status"]="complete";
  metadata["reference_present"]=normal_counts[o.ge_reference]>0;
  metadata["channel_coverage"]=coverage;
  metadata["histograms"]=writer.descriptions;
  metadata["source"]={{"path",o.input.string()},{"bytes",before.value.st_size},{"inode",before.value.st_ino},
    {"mtime_sec",before.value.st_mtim.tv_sec},{"mtime_nsec",before.value.st_mtim.tv_nsec},
    {"tree","tree"},{"total_entries",total},{"entry_begin",o.begin},{"entry_end",end},
    {"full_run",o.begin==0 && end==total},{"read_fields_sha256",digest.finish()},
    {"hash_definition","sr,pileup,outofr,cid,sid,ch,evte,ts,cfd,cfdft,cfds in entry order; native little-endian scalar bytes; no padding"}};
  metadata["counters"]={{"read",count},{"rejected_quality",rejected_quality},{"rejected_energy",rejected_energy},
      {"unmapped",unmapped},{"peak_buffer_hits",engine.peak_buffer}};
  metadata["chunk_entries"]=o.chunk;
  metadata["seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  writer.write(); output->cd();
  TObjString record(metadata.dump().c_str()); require(record.Write("metadata_json")>0,"Cannot write metadata");
  output->Flush(); require(!output->TestBit(TFile::kWriteError),"ROOT write error"); output->Close();
  { std::unique_ptr<TFile> check(TFile::Open(partial.c_str(),"READ"));
    require(check && !check->IsZombie(),"Cannot read back output");
    auto* m=dynamic_cast<TObjString*>(check->Get("metadata_json"));
    require(m && json::parse(m->GetString().Data()).at("status")=="complete","Output metadata readback failed");
    for(const auto& kv:writer.hist) {
      auto* h=dynamic_cast<TH1D*>(check->Get(kv.first.c_str()));
      require(h && h->GetEntries()==kv.second.all->GetEntries(),"Histogram readback failed");
    }
  }
  // link 只在目标不存在时成功，防止核验过程中出现的新文件被覆盖。
  fs::create_hard_link(partial,o.output); fs::remove(partial);
  std::cout<<json({{"status","complete"},{"output",o.output.string()},{"read_entries",count},
    {"reference_present",metadata["reference_present"]},{"seconds",metadata["seconds"]}}).dump()<<'\n';
}
