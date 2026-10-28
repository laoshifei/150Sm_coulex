// Inspect stored build timestamps without applying calibration or timing again.
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TObjString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TUUID.h>
#include <nlohmann/json.hpp>
#include <openssl/evp.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using json = nlohmann::json;
namespace fs = std::filesystem;
static void require(bool ok, const std::string &message) {
  if (!ok) throw std::runtime_error(message);
}
static std::string sha256(const fs::path &path) {
  std::ifstream in(path, std::ios::binary);
  require(bool(in), "Cannot hash " + path.string());
  std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> ctx(EVP_MD_CTX_new(), EVP_MD_CTX_free);
  require(ctx && EVP_DigestInit_ex(ctx.get(), EVP_sha256(), nullptr) == 1, "SHA initialization failed");
  std::vector<char> buffer(8 * 1024 * 1024);
  while (in) {
    in.read(buffer.data(), buffer.size());
    require(EVP_DigestUpdate(ctx.get(), buffer.data(), in.gcount()) == 1, "SHA update failed");
  }
  require(in.eof(), "Input error while hashing " + path.string());
  unsigned char bytes[EVP_MAX_MD_SIZE]; unsigned int size = 0;
  require(EVP_DigestFinal_ex(ctx.get(), bytes, &size) == 1, "SHA finalization failed");
  std::ostringstream out;
  for (unsigned int i = 0; i < size; ++i) out << std::hex << std::setw(2) << std::setfill('0') << int(bytes[i]);
  return out.str();
}
static json object(TFile &file, const char *name) {
  auto *value = dynamic_cast<TObjString *>(file.Get(name));
  require(value, std::string("Missing JSON object: ") + name);
  return json::parse(value->GetString().Data());
}
struct Family {
  TTreeReaderArray<Short_t> sid, ch, side;
  TTreeReaderArray<Double_t> energy, dt, time, alpha;
  TTreeReaderArray<Long64_t> source;
  Family(TTreeReader &r, const std::string &n)
      : sid(r,(n+"_sid").c_str()), ch(r,(n+"_ch").c_str()), side(r,(n+"_random_side").c_str()),
        energy(r,(n+"_energy").c_str()), dt(r,(n+"_dt_ns").c_str()), time(r,(n+"_time_ns").c_str()),
        alpha(r,(n+"_alpha").c_str()), source(r,(n+"_source_entry").c_str()) {}
  size_t size() const { return dt.GetSize(); }
  int key(size_t i) const { return 10000+100*sid[i]+ch[i]; }
  void validate() const {
    const size_t n=size();
    require(sid.GetSize()==n && ch.GetSize()==n && side.GetSize()==n && energy.GetSize()==n &&
      time.GetSize()==n && alpha.GetSize()==n && source.GetSize()==n, "Inconsistent particle arrays");
  }
};
struct Histograms {
  std::map<std::string,std::unique_ptr<TH1D>> one;
  std::map<std::string,std::unique_ptr<TH2D>> two;
  TH1D &h1(const std::string &name) {
    auto &p=one[name];
    if (!p) { p=std::make_unique<TH1D>(name.c_str(),";Corrected time difference (ns);Pairs / 2 ns",1001,-1001,1001); p->SetDirectory(nullptr); }
    return *p;
  }
  TH2D &h2(const std::string &name) {
    auto &p=two[name];
    if (!p) { p=std::make_unique<TH2D>(name.c_str(),";Ge energy (keV);Particle - Ge (ns)",98,40,2000,501,-1002,1002); p->SetDirectory(nullptr); }
    return *p;
  }
  void pair(const std::string &base, double dt, double energy, bool failed, int side) {
    const std::string region=side==0 ? "prompt" : (side<0 ? "random_left" : "random_right");
    const std::string state=failed ? "failed" : "normal";
    h1(base+"_"+region).Fill(dt);
    h1(base+"_"+region+"_"+state).Fill(dt);
    if (!side) h2(base+"_energy_"+state).Fill(energy,dt);
  }
  void write(TFile &out) {
    out.cd();
    for (auto &p:one) require(p.second->Write()>0,"Histogram write failed");
    for (auto &p:two) require(p.second->Write()>0,"Histogram write failed");
  }
};
struct Gate {
  double low, high, alpha;
  std::vector<std::pair<double,double>> random;
  bool accepts(double dt,int side) const {
    if (side==0) return low<=dt && dt<high;
    for (const auto &range:random)
      if ((range.second<=low ? -1 : 1)==side && range.first<=dt && dt<range.second) return true;
    return false;
  }
};
static std::map<std::pair<int,int>,Gate> gates(const json &manifest) {
  std::map<std::pair<int,int>,Gate> result;
  for (const auto &row:manifest.at("windows")) {
    Gate g{row.at("lower_ns"),row.at("upper_ns"),row.at("alpha"),{}};
    double width=0;
    for (const auto &range:row.at("random")) {
      const double a=range.at(0),b=range.at(1);
      require(a<b && (b<=g.low || a>=g.high),"Invalid embedded random interval");
      g.random.emplace_back(a,b); width+=b-a;
    }
    require(g.low<g.high && width>0 && std::abs(g.alpha-(g.high-g.low)/width)<1e-12,"Invalid embedded alpha");
    require(result.emplace(std::make_pair(row.at("ge").get<int>(),row.at("reference").get<int>()),g).second,"Duplicate embedded gate");
  }
  return result;
}
static double peak(const std::array<Long64_t,1001> &counts) {
  int best=2; double height=-1;
  for (int i=2;i<=998;++i) {
    double value=0;
    for(int j=-2;j<=2;++j) value+=counts[i+j]*(3-std::abs(j));
    if(value>height) { height=value;best=i; }
  }
  if(height<=0) return std::numeric_limits<double>::quiet_NaN();
  return -1000+2*best;
}
static json inspect(const fs::path &path, Histograms &hist, json &common) {
  const auto size=fs::file_size(path);
  const auto mtime=fs::last_write_time(path);
  const std::string input_sha=sha256(path);
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(),"READ"));
  require(file && !file->IsZombie(),"Cannot open "+path.string());
  const auto m=object(*file,"sm150_build_manifest"),expected=object(*file,"sm150_build_counts");
  require(m.at("schema_version")==2,"Only current schema-2 build files are supported");
  require(m.at("selection").at("single_Ge")==true,"Input is not the single-Ge production build");
  const int run=std::stoi(m.at("run").get<std::string>());
  const auto receipt_path=path.parent_path()/("data"+m.at("run").get<std::string>()+"_build_cpp_complete.json");
  json receipt; std::ifstream receipt_in(receipt_path); require(bool(receipt_in),"Missing production receipt"); receipt_in>>receipt;
  bool hash_ok=false;
  for(const auto &o:receipt.at("outputs"))
    if(fs::path(o.at("path").get<std::string>()).filename()==path.filename())
      hash_ok=o.at("sha256")==input_sha && o.at("bytes")==size;
  require(hash_ok,"ROOT SHA/size disagrees with production receipt: "+path.string());
  json identity;
  for(const auto &key:{"group","candidate","selection","time_alignment_profile","time_alignment_calibration_sha256","windows"}) identity[key]=m.at(key);
  if(common.is_null()) common=identity;
  require(common==identity,"Do not mix beam periods, candidates, timing versions or gate definitions");
  auto gate=gates(m);
  for(const auto &item:gate) {
    const auto base="ge"+std::to_string(item.first.first)+(item.first.second==0 ? "_si" : "_s3");
    for(const auto &region:{"prompt","random_left","random_right"}) {
      hist.h1(base+"_"+region);
      hist.h1(base+"_"+region+"_normal"); hist.h1(base+"_"+region+"_failed");
    }
    hist.h2(base+"_energy_normal");hist.h2(base+"_energy_failed");
  }
  // Source-entry identity removes repeat particle pairs associated with different Ge hits.
  std::set<std::tuple<int,Long64_t,Long64_t>> particle_pairs;
  std::map<std::string,std::array<Long64_t,1001>> per_run;
  json stats={{"run",run},{"path",fs::absolute(path).string()},{"sha256",input_sha},{"bytes",size},
    {"input_uuid",file->GetUUID().AsString()},{"particle_pair_repeats_removed",0},{"channel_mapping",m.at("energy_parameters")}};
  double max_dt_error=0;
  for(const auto &name:{"tr_event","tr_bg"}) {
    auto *tree=dynamic_cast<TTree *>(file->Get(name)); require(tree,"Missing build tree");
    const bool prompt=std::string(name)=="tr_event";
    for(const auto &branch:{"run","n_ge","ge_sid","ge_ch","ge_energy","ge_time_ns","ge_flags","ge_source_entry"})
      require(tree->GetBranch(branch),std::string("Missing branch ")+branch);
    for(const auto &family:{"spider","s3_ring","s3_sector"})
      for(const auto &field:{"sid","ch","random_side","energy","dt_ns","time_ns","alpha","source_entry"})
        require(tree->GetBranch((std::string(family)+"_"+field).c_str()),"Missing particle branch");
    TTreeReader reader(tree);
    TTreeReaderValue<Int_t> input_run(reader,"run"),n_ge(reader,"n_ge");
    TTreeReaderArray<Short_t> ge_sid(reader,"ge_sid"),ge_ch(reader,"ge_ch");
    TTreeReaderArray<Double_t> ge_energy(reader,"ge_energy"),ge_time(reader,"ge_time_ns");
    TTreeReaderArray<UChar_t> ge_flags(reader,"ge_flags");
    TTreeReaderArray<Long64_t> ge_source(reader,"ge_source_entry");
    Family si(reader,"spider"),ring(reader,"s3_ring"),sector(reader,"s3_sector");
    Long64_t rows=0,si_count=0,s3_count=0;
    while(reader.Next()) {
      ++rows;
      require(*input_run==run && *n_ge==1 && ge_sid.GetSize()==1 && ge_ch.GetSize()==1 && ge_energy.GetSize()==1 &&
        ge_time.GetSize()==1 && ge_flags.GetSize()==1 && ge_source.GetSize()==1,"Invalid Ge event");
      const int ge=100*ge_sid[0]+ge_ch[0]; const bool failed=(ge_flags[0]&4)!=0;
      require(std::isfinite(ge_time[0]) && ge_energy[0]>=40 && ge_energy[0]<2000 && !(ge_flags[0]&3),"Invalid Ge selection");
      si.validate();ring.validate();sector.validate();
      require(ring.size()==sector.size() && ring.size()<=1,"S3 ring/sector matching mismatch");
      for(int view=0;view<2;++view) {
        auto &family=view==0 ? si : ring;
        const auto base="ge"+std::to_string(ge)+(view==0 ? "_si" : "_s3");
        for(size_t i=0;i<family.size();++i) {
          const double dt=family.dt[i]; const int side=family.side[i];
          require(std::isfinite(dt) && std::isfinite(family.time[i]),"Nonfinite timestamp");
          require((prompt && side==0) || (!prompt && (side==-1 || side==1)),"Incorrect random-side flag");
          const auto &g=gate.at({ge,view});
          require(g.accepts(dt,side) && std::abs(g.alpha-family.alpha[i])<1e-12,"Pair outside embedded gate/alpha");
          require(family.energy[i]>(view==0 ? 15000 : 25000) && family.energy[i]<100000,"Particle outside stored energy gate");
          const double error=std::abs((family.time[i]-ge_time[0])-dt);
          max_dt_error=std::max(max_dt_error,error);
          const double tolerance=std::max(1e-6,4*std::numeric_limits<double>::epsilon()*std::abs(ge_time[0]));
          require(error<=tolerance,"Stored dt is inconsistent with corrected timestamps");
          hist.pair(base,dt,ge_energy[0],failed,side);
          if(prompt) {
            hist.h1("channel"+std::to_string(family.key(i))+"_ge_prompt").Fill(dt);
            const int bin=int(std::floor((dt+1001)/2));
            if(bin>=0 && bin<1001) { ++per_run[base][bin]; ++per_run[base+(failed ? "_failed" : "_normal")][bin]; }
          }
          view==0 ? ++si_count : ++s3_count;
        }
      }
      if(!prompt) continue;
      auto pair=[&](int kind,Long64_t a,Long64_t b,const std::string &label,double dt) {
        require(a!=b,"Self-pair in particle control");
        if(!particle_pairs.emplace(kind,a,b).second) {
          stats["particle_pair_repeats_removed"]=stats["particle_pair_repeats_removed"].get<Long64_t>()+1; return false;
        }
        hist.h1(label).Fill(dt); return true;
      };
      if(ring.size()) {
        const double dt=ring.dt[0]-sector.dt[0];
        require(std::abs(dt)<80 && std::abs(ring.energy[0]-sector.energy[0])<2500,"S3 match violates production cuts");
        if(pair(0,ring.source[0],sector.source[0],"s3_ring_minus_sector",dt)) {
          hist.h1("channel"+std::to_string(ring.key(0))+"_ring_minus_sector").Fill(dt);
          hist.h1("channel"+std::to_string(sector.key(0))+"_ring_minus_sector").Fill(dt);
        }
        for(size_t i=0;i<si.size();++i)
          if(pair(1,si.source[i],ring.source[0],"si_minus_s3",si.dt[i]-ring.dt[0]))
            hist.h1("channel"+std::to_string(si.key(i))+"_si_minus_s3").Fill(si.dt[i]-ring.dt[0]);
      }
      for(size_t i=0;i<si.size();++i) for(size_t j=i+1;j<si.size();++j) {
        if(si.key(i)==si.key(j)) continue;
        const size_t a=si.key(i)<si.key(j) ? i : j,b=a==i ? j : i;
        pair(2,si.source[a],si.source[b],"si_low_channel_minus_high_channel",si.dt[a]-si.dt[b]);
      }
    }
    require(rows==tree->GetEntries() && reader.GetEntryStatus()==TTreeReader::kEntryBeyondEnd,"Incomplete ROOT read");
    require(rows==expected.at(name).at("events") && si_count==expected.at(name).at("si_pairs") && s3_count==expected.at(name).at("s3_pairs"),"Read counts disagree with build metadata");
    stats[name]={{"events",rows},{"si_pairs",si_count},{"s3_pairs",s3_count}};
  }
  for(const auto &name:{"s3_ring_minus_sector","si_minus_s3","si_low_channel_minus_high_channel"}) hist.h1(name);
  stats["max_timestamp_subtraction_error_ns"]=max_dt_error;
  stats["per_ge_prompt"]=json::object();
  for(const auto &p:per_run) {
    Long64_t count=0;for(const auto n:p.second) count+=n;
    const double center=peak(p.second);
    stats["per_ge_prompt"][p.first]={{"pairs",count},{"local_mode_ns",std::isfinite(center) ? json(center) : json(nullptr)}};
  }
  require(fs::file_size(path)==size && fs::last_write_time(path)==mtime,"Input changed during inspection");
  std::cout<<"run "<<run<<" "<<m.at("group")<<" "<<m.at("candidate")<<" rows "<<stats["tr_event"]["events"]<<" + "<<stats["tr_bg"]["events"]<<std::endl;
  return stats;
}
int main(int argc,char **argv) {
  try {
    fs::path list,output;
    for(int i=1;i<argc;++i) {
      const std::string arg=argv[i];
      if(arg=="--help") { std::cout<<"check_build --input-list files.txt --output histograms.root\nAll input files must share a period, candidate and parameter version.\n";return 0; }
      require(i+1<argc,"Missing CLI argument");
      if(arg=="--input-list") list=argv[++i];else if(arg=="--output") output=argv[++i];else throw std::runtime_error("Unknown option: "+arg);
    }
    require(!list.empty() && !output.empty(),"Use --input-list and --output");
    const fs::path partial=output.string()+".partial";
    require(!fs::exists(output) && !fs::exists(partial),"Refuse to overwrite output/partial");
    std::ifstream in(list);require(bool(in),"Cannot open input list");
    std::vector<fs::path> paths;std::set<fs::path> seen;std::string line;
    while(std::getline(in,line)) {
      if(line.empty() || line[0]=='#') continue;
      if(line.back()=='\r') line.pop_back();
      const auto path=fs::canonical(line);
      require(seen.insert(path).second,"Duplicate input file");paths.push_back(path);
    }
    require(!paths.empty(),"Empty input list");
    TH1::AddDirectory(false);
    Histograms hist;json common,runs=json::array();
    for(const auto &path:paths) runs.push_back(inspect(path,hist,common));
    if(!output.parent_path().empty()) fs::create_directories(output.parent_path());
    std::unique_ptr<TFile> out(TFile::Open(partial.c_str(),"CREATE"));require(out && !out->IsZombie(),"Cannot create output");
    out->SetCompressionSettings(504);hist.write(*out);
    json metadata={{"schema_version",1},{"program","check_ts_align/check_build"},{"common",common},{"runs",runs},
      {"time_bin_ns",2},{"energy_time_bins",{{"energy_keV",20},{"time_ns",4}}},
      {"time_corrections_applied",false},{"prompt_and_random_are_disjoint_accepted_regions",true},
      {"particle_controls","unique source pairs within prompt rows only; Si-Si sign is lower hardware channel minus higher hardware channel"},
      {"limitations","Gated build events cannot recover full tails or independent ungated particle controls. S3 ring-sector already selected at abs(dt)<80 ns."}};
    require(TObjString(metadata.dump().c_str()).Write("sm150_timing_check_manifest")>0,"Metadata write failed");
    // Histograms were written explicitly and detached, so a final Write() may return zero.
    require(out->Write()>=0,"Output write failed");out->Close();require(!out->TestBit(TFile::kWriteError),"Output close failed");
    {
      std::unique_ptr<TFile> check(TFile::Open(partial.c_str(),"READ"));
      require(check && !check->IsZombie() && object(*check,"sm150_timing_check_manifest")==metadata,"Output metadata readback failed");
      for(const auto &p:hist.one) {
        auto *saved=dynamic_cast<TH1D *>(check->Get(p.first.c_str()));
        require(saved && saved->GetEntries()==p.second->GetEntries(),"Output histogram readback failed");
      }
      for(const auto &p:hist.two) {
        auto *saved=dynamic_cast<TH2D *>(check->Get(p.first.c_str()));
        require(saved && saved->GetEntries()==p.second->GetEntries(),"Output energy histogram readback failed");
      }
    }
    fs::create_hard_link(partial,output);fs::remove(partial);
    std::cout<<"wrote "<<output<<" from "<<paths.size()<<" complete input files\n";
    return 0;
  } catch(const std::exception &e) { std::cerr<<"check_build: "<<e.what()<<'\n';return 1; }
}
