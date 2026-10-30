#include "analysis.h"
#include "s3_tree.h"
#include "sm150_event_calibration.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TObjString.h>
#include <TTree.h>
#include <TUUID.h>
#include <chrono>
#include <iostream>
#include <memory>
#include <set>
#include <unistd.h>
using namespace sm150::diagnostic;

namespace {
constexpr double pi=3.14159265358979323846;
struct Pair { S3Hit ring,sector; double energy() const{return (ring.energy+sector.energy)/2.;} };
struct Geometry {
  std::vector<double> theta=std::vector<double>(25),phi=std::vector<double>(33);
  json snapshots=json::object();
  void read(const fs::path& path,std::vector<double>& table,double maximum) {
    auto text=read_text(path);snapshots[path.string()]={{"sha256",hash(text)},{"text",text}};
    std::istringstream lines(text);std::string line;std::getline(lines,line);std::set<int> seen;
    while(std::getline(lines,line)) {
      const auto comment=line.find('#');if(comment!=std::string::npos)line.resize(comment);
      std::istringstream row(line);int id;double angle;std::string extra;
      if(!(row>>id)){require(row.eof(),"Invalid geometry row");continue;}
      require(bool(row>>angle)&&!(row>>extra),"Invalid geometry row");
      require(id>0&&size_t(id)<table.size()&&seen.insert(id).second&&std::isfinite(angle)&&angle>=0&&angle<=maximum,"Invalid geometry value");
      table[id]=angle;
    }
    require(seen.size()+1==table.size(),"Incomplete geometry table");
  }
};
struct Energy {
  bool enabled=false;json snapshot;
  std::vector<double> a,b,q;
  std::map<int,Channel> base;
  Energy(const AnalysisOptions& o,const json& metadata) {
    if(o.energy_file.empty())return;
    enabled=true;const auto text=read_text(o.energy_file);auto all=json::parse(text);
    require(all.at("schema")=="sm150_s3_targeted_review_candidate_v1","Unsupported candidate schema");
    auto row=all.at("runs").at(rn(o.run));
    require(row.at("period")==metadata.at("period"),"Candidate period mismatch");
    a=row.at("A_keV").get<std::vector<double>>();b=row.at("B").get<std::vector<double>>();q=row.at("q_position").get<std::vector<double>>();
    require(a.size()==56&&b.size()==56&&q.size()==24,"Incomplete A/B/q candidate");
    for(size_t i=0;i<56;++i)require(std::isfinite(a[i])&&std::isfinite(b[i])&&b[i]>0,"Invalid A/B coefficient");
    for(double x:q)require(std::isfinite(x)&&x>0,"Invalid position coefficient");
    for(auto it=metadata.at("channels").begin();it!=metadata.at("channels").end();++it) {
      const auto& v=it.value();int kind=v.at("kind");if(kind!=3&&kind!=4)continue;
      Channel c;c.kind=kind;c.ring=v.at("ring");c.sector=v.at("sector");
      c.cal={v.at("cal_p0_kev"),v.at("cal_p1_kev_per_adc"),v.at("cal_p2_kev_per_adc2")};
      c.gain=v.at("gain");require(v.at("k_ge")==1.0,"Unexpected S3 Ge gain");base.emplace(std::stoi(it.key()),c);
    }
    require(base.size()==56,"Incomplete ADC calibration snapshot");
    snapshot={{"path",fs::absolute(o.energy_file).string()},{"sha256",hash(text)},{"run_parameters",row}};
  }
  double value(const S3Hit& hit,bool sector,int ring) const {
    if(!enabled)return hit.energy;
    const auto& c=base.at(hit.key);
    require(c.kind==(sector?4:3)&&(sector?c.sector:c.ring)==hit.id,"ADC mapping mismatch");
    const double e0=c.energy(hit.adc,false);
    const int index=(sector?24:0)+hit.id-1;
    return q.at(ring-1)*(a.at(index)+b.at(index)*e0);
  }
  Pair pair(const S3Hit& ring,const S3Hit& sector) const {
    Pair p{ring,sector};p.ring.energy=value(ring,false,ring.id);p.sector.energy=value(sector,true,ring.id);
    require(std::isfinite(p.ring.energy)&&std::isfinite(p.sector.energy),"Nonfinite corrected energy");return p;
  }
};
struct Histograms {
  std::map<std::string,std::unique_ptr<TH1>> hist;
  void one(const char* n,const char* title,int bins,double lo,double hi) {
    auto h=std::make_unique<TH1D>(n,title,bins,lo,hi);h->SetDirectory(nullptr);hist.emplace(n,std::move(h));
  }
  void two(const char* n,const char* title,int nx,double x0,double x1,int ny,double y0,double y1) {
    auto h=std::make_unique<TH2D>(n,title,nx,x0,x1,ny,y0,y1);h->SetDirectory(nullptr);hist.emplace(n,std::move(h));
  }
  void fill(const std::string& n,double x){hist.at(n)->Fill(x);}
  void fill(const std::string& n,double x,double y){static_cast<TH2*>(hist.at(n).get())->Fill(x,y);}
};
bool adjacent(int a,int b,int wrap=0) {return std::abs(a-b)==1||(wrap&&std::abs(a-b)==wrap-1);}
bool absent(const fs::path& p){return !fs::exists(p)&&!fs::is_symlink(p);}
}

int analysis::Process() {
  const auto started=std::chrono::steady_clock::now();
  auto o=options_;
  require(o.run>0&&o.run<=9999&&!o.input.empty()&&!o.output.empty(),"Specify RUN, --input and --output; see --help");
  require(o.segments>0&&o.segments<=100000&&o.sample_entries>=0,"Invalid sampling options");
  for(double v:{o.energy_min,o.energy_max,o.face_time,o.face_energy,o.particle_time,o.phi_low,o.phi_high,o.cm_low,o.cm_high,o.beam,o.deficit,o.projectile_mass,o.target_mass})
    require(std::isfinite(v),"Nonfinite option");
  require(o.energy_min>=0&&o.energy_max>o.energy_min&&o.face_time>0&&o.face_energy>0&&o.particle_time>0,"Invalid energy/time gate");
  require(o.phi_low>=0&&o.phi_high>o.phi_low&&o.phi_high<=360&&o.cm_low>=0&&o.cm_high>o.cm_low&&o.cm_high<=360,"Invalid angle gate");
  require(o.beam>0&&o.deficit>=0&&o.projectile_mass>0&&o.target_mass>0,"Invalid kinematic hypothesis");
  const double effective=o.beam-o.deficit*(1+o.projectile_mass/o.target_mass);
  require(effective>0,"Nonpositive effective beam energy");
  const double tau_r=std::sqrt(o.beam/effective),tau_p=o.projectile_mass/o.target_mass*tau_r;
  if(o.project.empty())o.project=fs::canonical("/proc/self/exe").parent_path().parent_path();
  o.input=fs::absolute(o.input);o.output=fs::absolute(o.output);o.project=fs::absolute(o.project);
  const auto partial=fs::path(o.output.string()+".partial"),receipt=fs::path(o.output.string()+".json");
  require(absent(o.output)&&absent(partial)&&absent(receipt),"Output or partial/receipt already exists; choose a new output path");
  const Stat source_before(o.input);
  std::unique_ptr<TFile> input(TFile::Open(o.input.c_str(),"READ"));
  require(input&&!input->IsZombie(),"Cannot open input ROOT");
  require(!input->Get("sm150_build_manifest"),"Ge-gated build_all ROOT retains at most one S3 pair; use build_s3 s3_events");
  require(o.format=="auto"||o.format=="s3_events"||o.format=="legacy","Invalid format");
  if(o.format=="auto") {
    require(input->Get("s3_events"),"No complete s3_events tree; legacy input requires explicit --format legacy");o.format="s3_events";
  }
  const bool modern=o.format=="s3_events";
  auto* tree=dynamic_cast<TTree*>(input->Get(modern?"s3_events":"tr_event"));
  require(tree&&!tree->GetBranch("n_ge"),"Missing S3 tree or unsupported Ge-gated input");
  const auto total=tree->GetEntries();require(total>0,"Empty input tree");
  json input_meta=json::object();std::string metadata_text;
  if(modern) {
    auto* object=dynamic_cast<TObjString*>(input->Get("metadata"));require(object,"Missing s3_events metadata");
    metadata_text=object->GetString().Data();input_meta=json::parse(metadata_text);
    require(input_meta.at("schema")=="sm150_s3_calibration_v1","Unsupported S3 schema");
    require(input_meta.at("run")==rn(o.run),"Input run mismatch");
  } else require(o.energy_file.empty(),"Legacy input lacks ADC snapshots for candidate energy");
  Geometry geometry;
  geometry.read(o.project/"pars/angle_info/s3_ring_angle.txt",geometry.theta,180.);
  geometry.read(o.project/"pars/angle_info/s3_sector_angle.txt",geometry.phi,360.);
  Energy energy(o,input_meta);
  const auto count=o.sample_entries==0?total:std::min<Long64_t>(total,o.sample_entries);
  const int segments=count==total?1:std::min<Long64_t>(o.segments,count);
  std::vector<std::pair<Long64_t,Long64_t>> ranges;
  json ranges_json=json::array();
  for(int i=0;i<segments;++i) {
    const auto begin=total/segments*i+(total%segments*i)/segments;
    const auto end=total/segments*(i+1)+(total%segments*(i+1))/segments;
    const auto length=count/segments+(i<count%segments?1:0),first=begin+(end-begin-length)/2;
    require(length>0&&first>=begin&&first+length<=end,"Invalid sample range");
    ranges.emplace_back(first,first+length);ranges_json.push_back({first,first+length});
  }
  s3_tree reader(tree,o.format);
  Histograms h;
  h.two("multiplicity_raw","Input clusters;ring multiplicity;sector multiplicity",16,-.5,15.5,16,-.5,15.5);
  h.two("multiplicity_clean","Unflagged interior clusters;ring multiplicity;sector multiplicity",16,-.5,15.5,16,-.5,15.5);
  h.two("ring_energy_single","Accepted 1+1;ring ID;ring energy [keV]",24,.5,24.5,200,0,100000);
  h.one("assignments","2+2 clusters;passing bijections;clusters",3,-.5,2.5);
  h.one("face_dt_all","All four 2+2 pair possibilities;ring - sector time [ns];candidate pairs",240,-300,300);
  h.one("face_de_all","All four 2+2 pair possibilities;ring - sector energy [keV];candidate pairs",400,-100000,100000);
  h.one("face_dt_matched","Unique 2+2 assignment;ring - sector time [ns];matched pairs",160,-80,80);
  h.one("face_de_matched","Unique 2+2 assignment;ring - sector energy [keV];matched pairs",200,-2500,2500);
  h.two("face_de_energy","Unique 2+2 assignment;mean pair energy [keV];ring - sector energy [keV]",100,0,100000,100,-2500,2500);
  h.one("particle_dt","Unique 2+2 assignment;high - low pair mean time [ns];clusters",240,-300,300);
  h.one("phi_before_time","Unique 2+2 assignment;absolute sector azimuth difference [deg];clusters",144,0,360);
  h.one("phi_after_time","After pair-time gate;absolute sector azimuth difference [deg];clusters",144,0,360);
  h.one("cm_before_phi","After pair-time gate;CM angle sum under legacy hypothesis [deg];clusters",180,0,360);
  h.one("cm_after_phi","After pair-time and azimuth gates;CM angle sum under legacy hypothesis [deg];clusters",180,0,360);
  for(const std::string tag:{"before","after"}) {
    h.two(("energy_low_high_"+tag).c_str(),"Two particles;low mean energy [keV];high mean energy [keV]",200,0,100000,200,0,100000);
    h.two(("theta_low_high_"+tag).c_str(),"Two particles;low-energy ring theta [deg];high-energy ring theta [deg]",90,0,90,90,0,90);
    h.two(("energy_theta_"+tag).c_str(),"Two particles;ring theta [deg];mean pair energy [keV]",90,0,90,200,0,100000);
    h.one(("energy_low_"+tag).c_str(),"Low-energy particle;mean pair energy [keV];clusters",1000,0,100000);
    h.one(("energy_high_"+tag).c_str(),"High-energy particle;mean pair energy [keV];clusters",1000,0,100000);
  }
  std::map<std::string,Long64_t> c;
  for(auto key:{"input","boundary","quality_rejected","clean","one_one","one_two","two_one","two_two","other","sharing_like","duplicate_channel","no_assignment","ambiguous","unique","particle_time","phi","cm_invalid","accepted","single_accepted","adjacent_ring","adjacent_sector"})c[key]=0;
  fs::create_directories(o.output.parent_path());
  std::unique_ptr<TFile> output(TFile::Open(partial.c_str(),"NEW"));
  require(output&&!output->IsZombie(),"Cannot create new output");output->SetCompressionSettings(504);
  TTree candidates("candidates","Unique 2+2 assignments, including those failing later angular/time gates");
  candidates.SetDirectory(output.get());candidates.SetAutoFlush(-8*1024*1024);
  Long64_t entry=0,cluster=0,rsrc[2]{},ssrc[2]{};Int_t run=o.run,rid[2]{},sid[2]{};
  Double_t er[2]{},es[2]{},tr[2]{},ts[2]{},dt=0,phi=0,cm=0;
  Bool_t adjacent_ring=false,adjacent_sector=false,cm_valid=false,pass_time=false,pass_phi=false,pass_cm=false;
  candidates.Branch("run",&run,"run/I");candidates.Branch("source_event_entry",&entry,"source_event_entry/L");candidates.Branch("cluster",&cluster,"cluster/L");
  candidates.Branch("ring_id",rid,"ring_id[2]/I");candidates.Branch("sector_id",sid,"sector_id[2]/I");
  candidates.Branch("ring_source_entry",rsrc,"ring_source_entry[2]/L");candidates.Branch("sector_source_entry",ssrc,"sector_source_entry[2]/L");
  candidates.Branch("ring_energy_keV",er,"ring_energy_keV[2]/D");candidates.Branch("sector_energy_keV",es,"sector_energy_keV[2]/D");
  candidates.Branch("ring_time_ns",tr,"ring_time_ns[2]/D");candidates.Branch("sector_time_ns",ts,"sector_time_ns[2]/D");
  candidates.Branch("particle_dt_ns",&dt,"particle_dt_ns/D");candidates.Branch("phi_difference_deg",&phi,"phi_difference_deg/D");
  candidates.Branch("cm_sum_deg",&cm,"cm_sum_deg/D");candidates.Branch("cm_valid",&cm_valid,"cm_valid/O");
  candidates.Branch("adjacent_ring",&adjacent_ring,"adjacent_ring/O");candidates.Branch("adjacent_sector",&adjacent_sector,"adjacent_sector/O");
  candidates.Branch("passed_particle_time",&pass_time,"passed_particle_time/O");candidates.Branch("passed_phi",&pass_phi,"passed_phi/O");candidates.Branch("passed_cm",&pass_cm,"passed_cm/O");
  auto valid=[&](const Pair& p) {return p.ring.energy>o.energy_min&&p.ring.energy<o.energy_max&&p.sector.energy>o.energy_min&&p.sector.energy<o.energy_max&&std::abs(p.ring.time-p.sector.time)<o.face_time&&std::abs(p.ring.energy-p.sector.energy)<o.face_energy;};
  Digest field_hash;
  for(const auto& range:ranges)for(Long64_t index=range.first;index<range.second;++index) {
    auto event=reader.Read(index);++c["input"];int nr=event.rings.size(),ns=event.sectors.size();
    field_hash.add(index);field_hash.add(event.cluster);field_hash.add(event.boundary);field_hash.add(nr);field_hash.add(ns);
    bool flagged=false;
    for(const auto* face:{&event.rings,&event.sectors})for(const auto& hit:*face) {
      flagged=flagged||hit.flags!=0;field_hash.add(hit.id);field_hash.add(hit.key);field_hash.add(hit.adc);field_hash.add(hit.flags);field_hash.add(hit.source);field_hash.add(hit.time);field_hash.add(hit.energy);
    }
    h.fill("multiplicity_raw",nr,ns);
    if(event.boundary){++c["boundary"];continue;}if(flagged){++c["quality_rejected"];continue;}
    ++c["clean"];h.fill("multiplicity_clean",nr,ns);
    if(nr==1&&ns==1) {
      ++c["one_one"];auto p=energy.pair(event.rings[0],event.sectors[0]);
      if(valid(p)){++c["single_accepted"];h.fill("ring_energy_single",p.ring.id,p.ring.energy);}continue;
    }
    if(nr==1&&ns==2){++c["one_two"];if(adjacent(event.sectors[0].id,event.sectors[1].id,32))++c["sharing_like"];continue;}
    if(nr==2&&ns==1){++c["two_one"];if(adjacent(event.rings[0].id,event.rings[1].id))++c["sharing_like"];continue;}
    if(nr!=2||ns!=2){++c["other"];continue;}++c["two_two"];
    if(event.rings[0].id==event.rings[1].id||event.sectors[0].id==event.sectors[1].id){++c["duplicate_channel"];continue;}
    std::array<std::array<Pair,2>,2> permutations;
    int passes=0,chosen=-1;
    for(int swap=0;swap<2;++swap) {
      bool good=true;
      for(int i=0;i<2;++i) {
        auto p=energy.pair(event.rings[i],event.sectors[i^swap]);permutations[swap][i]=p;
        h.fill("face_dt_all",p.ring.time-p.sector.time);h.fill("face_de_all",p.ring.energy-p.sector.energy);good=good&&valid(p);
      }
      if(good){++passes;chosen=swap;}
    }
    h.fill("assignments",passes);
    if(passes==0){++c["no_assignment"];continue;}if(passes==2){++c["ambiguous"];continue;}++c["unique"];
    auto pair=permutations[chosen];
    std::sort(pair.begin(),pair.end(),[](const Pair& a,const Pair& b){return a.energy()==b.energy()?a.ring.id<b.ring.id:a.energy()<b.energy();});
    entry=event.entry;cluster=event.cluster;
    for(int i=0;i<2;++i) {
      const auto& p=pair[i];rid[i]=p.ring.id;sid[i]=p.sector.id;rsrc[i]=p.ring.source;ssrc[i]=p.sector.source;
      er[i]=p.ring.energy;es[i]=p.sector.energy;tr[i]=p.ring.time;ts[i]=p.sector.time;
      h.fill("face_dt_matched",tr[i]-ts[i]);h.fill("face_de_matched",er[i]-es[i]);h.fill("face_de_energy",p.energy(),er[i]-es[i]);
    }
    adjacent_ring=adjacent(rid[0],rid[1]);adjacent_sector=adjacent(sid[0],sid[1],32);
    c["adjacent_ring"]+=adjacent_ring;c["adjacent_sector"]+=adjacent_sector;
    dt=((tr[1]-tr[0])+(ts[1]-ts[0]))/2.;phi=std::abs(geometry.phi[sid[1]]-geometry.phi[sid[0]]);
    h.fill("particle_dt",dt);h.fill("phi_before_time",phi);
    pass_time=std::abs(dt)<o.particle_time;pass_phi=phi>o.phi_low&&phi<o.phi_high;
    const double t0=geometry.theta[rid[0]],t1=geometry.theta[rid[1]],a0=tau_r*std::sin(t0*pi/180.),a1=tau_p*std::sin(t1*pi/180.);
    cm_valid=std::abs(a0)<=1&&std::abs(a1)<=1;cm=0;
    if(cm_valid)cm=t0+t1+(std::asin(a0)+std::asin(a1))*180./pi;else ++c["cm_invalid"];
    pass_cm=cm_valid&&cm>o.cm_low&&cm<o.cm_high;
    require(candidates.Fill()>0,"Candidate tree write failed");
    if(!pass_time)continue;
    ++c["particle_time"];h.fill("phi_after_time",phi);
    if(cm_valid)h.fill("cm_before_phi",cm);
    auto fill_pair=[&](const std::string& tag) {
      h.fill("energy_low_high_"+tag,pair[0].energy(),pair[1].energy());
      h.fill("theta_low_high_"+tag,t0,t1);h.fill("energy_low_"+tag,pair[0].energy());h.fill("energy_high_"+tag,pair[1].energy());
      h.fill("energy_theta_"+tag,t0,pair[0].energy());h.fill("energy_theta_"+tag,t1,pair[1].energy());
    };
    fill_pair("before");
    if(!pass_phi)continue;
    ++c["phi"];if(cm_valid)h.fill("cm_after_phi",cm);
    if(!pass_cm)continue;
    ++c["accepted"];fill_pair("after");
  }
  require(c["input"]==count&&candidates.GetEntries()==c["unique"],"Read/count mismatch");
  require(Stat(o.input).same(source_before),"Input file changed during analysis");
  for(auto it=geometry.snapshots.begin();it!=geometry.snapshots.end();++it)
    require(hash(read_text(it.key()))==it.value().at("sha256"),"Geometry changed during analysis");
  if(energy.enabled)require(hash(read_text(o.energy_file))==energy.snapshot.at("sha256"),"Energy candidate changed");
  const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
  json metadata={{"schema","sm150_s3_two_particle_v1"},{"run",rn(o.run)},{"format",o.format},
    {"input",o.input.string()},{"input_bytes",source_before.value.st_size},{"input_uuid",input->GetUUID().AsString()},
    {"input_entries",total},{"entries_read",count},{"sample_fraction",double(count)/total},{"segments",ranges_json},
    {"read_fields_sha256",field_hash.finish()},{"input_metadata",input_meta},{"input_metadata_sha256",hash(metadata_text)},
    {"energy_mode",energy.enabled?"explicit_complete_candidate":"stored_previous_energy"},
    {"energy_candidate",energy.snapshot},{"time_mode","stored corrected time; no extra correction"},
    {"geometry",geometry.snapshots},{"counts",c},{"elapsed_seconds",elapsed},
    {"gates",{{"energy_min_keV",o.energy_min},{"energy_max_keV",o.energy_max},{"face_time_ns",o.face_time},{"face_energy_keV",o.face_energy},{"particle_time_ns",o.particle_time},{"phi_low_deg",o.phi_low},{"phi_high_deg",o.phi_high},{"cm_low_deg",o.cm_low},{"cm_high_deg",o.cm_high},{"bounds","strict"}}},
    {"kinematic_hypothesis",{{"beam_MeV",o.beam},{"deficit_MeV",o.deficit},{"projectile_A",o.projectile_mass},{"target_A",o.target_mass},{"tau_recoil",tau_r},{"tau_projectile",tau_p},{"assignment","lower mean energy = recoil, higher = projectile; diagnostic hypothesis, not particle identification"}}},
    {"limitations",{"Only input two-sided connected time clusters are observable; no claim about lost single-face hits or wide random sidebands.","Adjacent strips are marked, not merged or vetoed; 1+2/2+1 adjacency alone does not establish charge sharing.","Both bijections are evaluated; two passing assignments are counted as ambiguous and excluded.","before histograms mean after unique assignment and pair-time gate, before angular gates.","No prompt-random subtraction; raw diagnostic counts only.","High-beam rings 9-24 are included for diagnostics, unlike production build selection."}}};
  output->cd();
  for(auto& kv:h.hist){require(kv.second->Write()>0,"Histogram write failed");metadata["histogram_entries"][kv.first]=kv.second->GetEntries();}
  require(candidates.Write()>0,"Tree write failed");
  TObjString object(metadata.dump().c_str());require(object.Write("metadata")>0,"Metadata write failed");
  candidates.SetDirectory(nullptr);output->Close();
  require(!output->TestBit(TFile::kWriteError),"ROOT reported a write error");
  {
    std::unique_ptr<TFile> check(TFile::Open(partial.c_str(),"READ"));require(check&&!check->IsZombie(),"Cannot reopen output");
    auto* saved=dynamic_cast<TTree*>(check->Get("candidates"));
    require(saved&&saved->GetEntries()==c["unique"]&&check->Get("metadata"),"Output readback failed");
    for(const auto& kv:h.hist){auto* hist=dynamic_cast<TH1*>(check->Get(kv.first.c_str()));require(hist&&hist->GetEntries()==kv.second->GetEntries(),"Histogram readback failed");}
  }
  require(::link(partial.c_str(),o.output.c_str())==0,"Cannot publish output without replacement");
  fs::remove(partial);
  std::ofstream summary(receipt);require(bool(summary),"Cannot create receipt");summary<<metadata.dump(2)<<'\n';summary.close();require(bool(summary),"Receipt write failed");
  std::cout<<metadata.at("counts").dump(2)<<"\nElapsed seconds: "<<elapsed<<"\nOutput: "<<o.output<<'\n';
  return 0;
}
