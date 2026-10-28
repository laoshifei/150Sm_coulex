#ifndef SM150_EVENT_CALIBRATION_H
#define SM150_EVENT_CALIBRATION_H
// Shared by get_ts_offset and ge_time_walk: retain validated calibration, mapping and parameter-identity checks.
// Convert ADC to double before evaluation; add no random time/energy jitter and modify no parameter files.
#include <nlohmann/json.hpp>
#include <openssl/evp.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>
namespace sm150 { namespace diagnostic {
namespace fs=std::filesystem;
using json=nlohmann::json;
inline void require(bool ok,const std::string& s) { if(!ok) throw std::runtime_error(s); }
inline std::string rn(int run) { std::ostringstream s; s<<std::setw(4)<<std::setfill('0')<<run; return s.str(); }
struct Digest {
  EVP_MD_CTX* context=EVP_MD_CTX_new();
  Digest() { require(context && EVP_DigestInit_ex(context,EVP_sha256(),nullptr)==1,"SHA256 init failed"); }
  ~Digest() { EVP_MD_CTX_free(context); }
  void add(const void* p,size_t n) { require(EVP_DigestUpdate(context,p,n)==1,"SHA256 failed"); }
  template<class T> void add(const T& x) { add(&x,sizeof(x)); }
  std::string finish() {
    unsigned char v[EVP_MAX_MD_SIZE]; unsigned n=0;
    require(EVP_DigestFinal_ex(context,v,&n)==1,"SHA256 final failed");
    std::ostringstream s; for(unsigned i=0;i<n;++i) s<<std::hex<<std::setw(2)<<std::setfill('0')<<unsigned(v[i]);
    return s.str();
  }
};
inline std::string hash(const std::string& text) { Digest d; d.add(text.data(),text.size()); return d.finish(); }
inline std::string read_text(const fs::path& p) {
  std::ifstream f(p,std::ios::binary); require(bool(f),"Cannot read "+p.string());
  std::string s((std::istreambuf_iterator<char>(f)),{});
  require(!f.bad(),"Read failed "+p.string()); return s;
}
struct Stat {
  struct stat value{};
  explicit Stat(const fs::path& p) { require(::stat(p.c_str(),&value)==0,"Cannot stat "+p.string()); }
  bool same(const Stat& b) const {
    return value.st_ino==b.value.st_ino && value.st_dev==b.value.st_dev &&
      value.st_size==b.value.st_size && value.st_mtim.tv_sec==b.value.st_mtim.tv_sec &&
      value.st_mtim.tv_nsec==b.value.st_mtim.tv_nsec;
  }
};
struct Channel {
  int kind=0, ring=0, sector=0; // 1 Ge, 2 Si, 3 S3 ring, 4 S3 sector
  std::array<double,3> cal{};
  double gain=1, kg=1, att=1, fb0=0, fb1=1;
  double energy(double adc,bool final) const {
    double e=(cal[0]+cal[1]*adc+cal[2]*adc*adc)*gain*kg*att;
    return final ? fb0+fb1*e : e;
  }
};
struct Parameters {
  fs::path pars;
  std::string group;
  std::map<int,Channel> channels;
  json snapshots=json::object(), records=json::object();
  std::vector<std::vector<double>> table(const fs::path& p,bool header=true) {
    const auto s=read_text(p); snapshots[p.string()]={{"sha256",hash(s)},{"text",s}};
    std::istringstream input(s); std::string line; if(header) std::getline(input,line);
    std::vector<std::vector<double>> rows;
    while(std::getline(input,line)) {
      auto pos=line.find('#'); if(pos!=std::string::npos) line.resize(pos);
      std::istringstream in(line); double x; std::vector<double> row;
      while(in>>x) { require(std::isfinite(x),"Nonfinite parameter "+p.string()); row.push_back(x); }
      require(in.eof(),"Non-numeric row "+p.string()); if(!row.empty()) rows.push_back(row);
    }
    require(!rows.empty(),"Empty table "+p.string()); return rows;
  }
  std::map<int,std::vector<double>> keyed(const fs::path& p,int count) {
    std::map<int,std::vector<double>> out;
    for(const auto& r:table(p)) {
      require(r.size()>=size_t(3+count),"Short parameter row "+p.string());
      require(r[0]==std::floor(r[0]) && r[1]==std::floor(r[1]) && r[2]==std::floor(r[2]) &&
        r[0]>=0 && r[0]<=1 && r[1]>=2 && r[1]<=12 && r[2]>=0 && r[2]<16,"Invalid channel row");
      int k=int(r[0])*10000+int(r[1])*100+int(r[2]);
      require(!out.count(k),"Duplicate parameter channel "+std::to_string(k));
      out[k]=std::vector<double>(r.begin()+3,r.begin()+3+count);
    }
    return out;
  }
  Parameters(const fs::path& project,int run,bool final):pars(project/"pars") {
    if(run>=472 && run<=605) group="high1";
    else if(run>=616 && run<=674) group="low";
    else if(run>=680 && run<=718) group="high2";
    else throw std::runtime_error("Only the three in-beam periods are supported");
    auto ge=keyed(pars/"cali_ge.txt",3);
    auto si=keyed(pars/(group=="low"?"cali_si_r457_low_current.txt":group=="high1"?"cali_si_r460.txt":"cali_si_r457_high_current.txt"),3);
    auto kg=keyed(pars/"run_k"/("k_"+rn(run)+".txt"),1);
    // Match the 37 physical Ge channels in build_all; do not implicitly remove channels based on config/b.txt.
    const int ge_keys[]={204,205,206,208,209,210,211,212,213,214,215,304,305,306,307,308,309,310,311,
                         400,401,402,403,404,405,406,407,409,410,411,504,505,506,507,508,509,510};
    auto add=[&](int key,int kind,int ring,int sector,double gain) {
      const auto& c=(kind==1?ge:si).at(key);
      Channel h; h.kind=kind; h.ring=ring; h.sector=sector; h.gain=gain;
      std::copy(c.begin(),c.begin()+3,h.cal.begin());
      if(kind==1) { h.kg=kg.at(key)[0]; require(h.kg>0,"Nonpositive Ge gain"); }
      require(channels.emplace(key,h).second,"Duplicate detector mapping");
    };
    for(int k:ge_keys) add(k,1,0,0,1);
    const int sa=group=="low"?9:2, sb=group=="low"?10:3;
    const int slots[]={sa,sb,4,6,7,8};
    for(int slot:slots) for(int ch=0;ch<16;++ch) {
      int ring=(slot==sa || slot==4)?8-ch%8:ch%8+1;
      int sector=slot==sa?(ch<8?2:1):slot==sb?(ch<8?4:3):slot==4?(ch<8?6:5):
                 slot==6?(ch<8?7:8):slot==7?(ch<8?9:10):(ch<8?11:12);
      double gain=((run>475 && run<616 && slot==2)||(run>=680 && slot==2)||(group=="low" && slot==9))?4.:1.;
      add(10000+slot*100+ch,2,ring,sector,gain);
    }
    for(int ch=2;ch<=13;++ch) add(11100+ch,3,2*(ch-2)+1,0,10.);
    for(int ch=0;ch<12;++ch) add(11200+ch,3,2*ch+2,0,1.);
    const int s0=group=="low"?2:9, s1=group=="low"?3:10;
    for(int ch=0;ch<16;++ch) {
      add(10000+100*s0+ch,4,0,(10+2*ch-1)%32+1,1.);
      add(10000+100*s1+ch,4,0,(11+2*ch-1)%32+1,1.);
    }
    for(int face=0;face<2;++face) {
      const std::string name=face?"sector":"ring"; const int count=face?32:24;
      std::map<int,std::array<double,2>> fb;
      if(final) for(const auto& r:table(pars/"run_fb_cor"/("correction_"+name+"_"+rn(run)+".txt"))) {
        require(r.size()>=3 && r[0]==std::floor(r[0]) && r[0]>=1 && r[0]<=count &&
                r[2]>0 && !fb.count(int(r[0])),"Invalid front-back table");
        fb[int(r[0])]={r[1],r[2]};
      }
      if(final) require(fb.size()==size_t(count),"Incomplete front-back table");
      for(auto& kv:channels) if(kv.second.kind==(face?4:3)) {
        auto& c=kv.second; int id=face?c.sector:c.ring;
        std::ostringstream number; number<<std::setw(2)<<std::setfill('0')<<id;
        int found=0;
        for(const auto& r:table(pars/"run_s3_att"/(name+number.str()+".txt"),false)) {
          require(r.size()>=2,"Short attenuation row");
          if(r[0]==run) { c.att=r[1]; ++found; }
        }
        require(found==1 && c.att>0,"Missing/invalid attenuation for "+name+number.str());
        if(final) { c.fb0=fb.at(id)[0]; c.fb1=fb.at(id)[1]; }
      }
    }
    require(channels.size()==189,"Incomplete physical mapping");
    for(const auto& kv:channels) {
      int key=kv.first; const auto& c=kv.second;
      records[std::to_string(key)]={{"kind",c.kind},{"ring",c.ring},{"sector",c.sector},
        {"cal_p0_kev",c.cal[0]},{"cal_p1_kev_per_adc",c.cal[1]},{"cal_p2_kev_per_adc2",c.cal[2]},
        {"gain",c.gain},{"k_ge",c.kg},{"attenuation",c.att},{"fb_p0_kev",c.fb0},{"fb_p1",c.fb1}};
    }
  }
  void freeze(const fs::path& p) { auto s=read_text(p); snapshots[p.string()]={{"sha256",hash(s)},{"text",s}}; }
  void verify() const {
    for(auto it=snapshots.begin();it!=snapshots.end();++it)
      require(hash(read_text(it.key()))==it.value().at("sha256"),"Parameter changed: "+it.key());
  }
};

} } // namespace sm150::diagnostic
#endif
