// Compact distributed S3 pairs for joint energy calibration.
// Every stored energy is computed before attenuation and front/back correction.
// Keep prompt, sidebands, multiplicity, and previous-hit information separately.
#include "build.h"
#include "sm150_event_calibration.h"
#include "sm150_time_alignment.h"
#include <TFile.h>
#include <TLeaf.h>
#include <TObjString.h>
#include <TTree.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>
#include <unistd.h>
namespace {
using namespace sm150::diagnostic;
struct Hit {
  Long64_t source=0,native=0,time=0,previous=-1;
  int channel=0,adc=0,flags=0;
  double energy=0,old=0,pre=0;
};
struct Group { Hit main; int count=0,flags=0; double minor=0; bool adjacent=false; };
std::vector<Group> group_face(const std::vector<Hit>& hits,Long64_t gap,int channels) {
  std::vector<Group> groups;
  for(size_t first=0;first<hits.size();) {
    size_t end=first+1,best=first;
    while(end<hits.size() && hits[end].time-hits[end-1].time<=gap)++end;
    for(size_t i=first+1;i<end;++i)if(hits[i].pre>hits[best].pre)best=i;
    Group g;g.main=hits[best];g.count=int(end-first);double other=0;
    for(size_t i=first;i<end;++i) {
      g.flags|=hits[i].flags;
      if(i==best)continue;
      other+=hits[i].pre;
      int separation=std::abs(hits[i].channel-hits[best].channel);
      if(separation==1 || (channels==32 && separation==31))g.adjacent=true;
    }
    g.minor=other/g.main.pre;groups.push_back(g);first=end;
  }
  std::sort(groups.begin(),groups.end(),[](const Group& a,const Group& b){return a.main.time<b.main.time;});
  return groups;
}
}
int build::Sample() {
  using namespace sm150::diagnostic;
  using Clock=std::chrono::steady_clock;
  try {
    const auto start=Clock::now();const auto& o=options_;
    require(o.blocks>=4 && o.blocks%2==0,"Joint sampling needs an even number of at least four blocks");
    const auto project=o.project.empty()?fs::canonical("/proc/self/exe").parent_path().parent_path():fs::canonical(o.project);
    const auto input=fs::canonical(o.input.empty()?fs::path("/mnt/e/Sm150/rootfile/data_C2_"+rn(o.run)+"_"+rn(o.run)+".root"):o.input);
    const auto output=fs::absolute(o.output.empty()?project/"rootfile"/("s3_joint_"+rn(o.run)+".root"):o.output);
    const fs::path partial=output.string()+".partial";
    require(!fs::exists(output) && !fs::exists(partial) && !fs::exists(output.string()+".json"),"Output or receipt already exists");
    Parameters pars(project,o.run,true);fs::path profile=o.profile_dir;
    if(profile.empty() && pars.paths.active()) {
      profile=pars.paths.timing_directory();
    } else if(profile.empty()) {
      const auto active=project/"pars/time_alignment/active.txt";pars.freeze(active);
      std::istringstream in(read_text(active));std::string name,extra;
      require(bool(in>>name) && !(in>>extra) && name.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_")==std::string::npos,"Invalid active timing profile");
      profile=project/"pars/time_alignment/versions"/name;
    }
    profile=fs::canonical(profile);
    for(const auto& f:fs::directory_iterator(profile))if(f.is_regular_file() && f.path().filename()!="README.md")pars.freeze(f.path());
    sm150::TimeAlignment timing;timing.LoadProfile(profile.string(),o.run);require(timing.total(),"Incomplete timing");
    pars.verify();const Stat before(input);
    std::unique_ptr<TFile> source(TFile::Open(input.c_str(),"READ"));require(source && !source->IsZombie(),"Cannot open C2");
    auto* tree=dynamic_cast<TTree*>(source->Get("tree"));require(tree,"No decoded tree");
    Long64_t total=tree->GetEntries(),wanted=o.sample_entries?std::min<Long64_t>(o.sample_entries,total):total;
    require(wanted>=o.blocks,"Too few sample entries");
    Short_t cid=0,sid=0,ch=0,sr=0;UShort_t adc=0;Long64_t ts=0;Bool_t pileup=false,outofr=false;
    const std::map<std::string,std::string> types={{"cid","Short_t"},{"sid","Short_t"},{"ch","Short_t"},{"sr","Short_t"},{"evte","UShort_t"},{"ts","Long64_t"},{"pileup","Bool_t"},{"outofr","Bool_t"}};
    tree->SetBranchStatus("*",false);
    for(const auto& kv:types) {
      auto* leaf=tree->GetLeaf(kv.first.c_str());require(leaf && kv.second==leaf->GetTypeName(),"Wrong branch: "+kv.first);
      tree->SetBranchStatus(kv.first.c_str(),true);
    }
    require(tree->SetBranchAddress("cid",&cid)>=0 && tree->SetBranchAddress("sid",&sid)>=0 && tree->SetBranchAddress("ch",&ch)>=0 &&
      tree->SetBranchAddress("sr",&sr)>=0 && tree->SetBranchAddress("evte",&adc)>=0 && tree->SetBranchAddress("ts",&ts)>=0 &&
      tree->SetBranchAddress("pileup",&pileup)>=0 && tree->SetBranchAddress("outofr",&outofr)>=0,"Cannot bind branches");
    tree->SetCacheSize(8*1024*1024);for(const auto& kv:types)tree->AddBranchToCache(kv.first.c_str(),true);tree->StopCacheLearningPhase();
    fs::create_directories(output.parent_path());
    std::unique_ptr<TFile> target(TFile::Open(partial.c_str(),"CREATE","S3 joint calibration sample",404));require(target && !target->IsZombie(),"Cannot create sample");
    TTree pairs("pairs","Primary face hits; no front/back energy-agreement selection");pairs.SetAutoFlush(-16*1024*1024);
    Int_t block=0,r=0,s=0,nr=0,ns=0,flags=0,group_flags=0,ra=0,sa=0;
    Bool_t adjacent_r=false,adjacent_s=false;
    Long64_t re=0,se=0,rt=0,st=0,rprev=0,sprev=0;
    Double_t er=0,es=0,oldr=0,olds=0,dt=0,mr=0,ms=0;
    pairs.Branch("block",&block);pairs.Branch("r",&r);pairs.Branch("s",&s);
    pairs.Branch("nr",&nr);pairs.Branch("ns",&ns);pairs.Branch("flags",&flags);pairs.Branch("group_flags",&group_flags);
    pairs.Branch("r_adc",&ra);pairs.Branch("s_adc",&sa);pairs.Branch("r_entry",&re);pairs.Branch("s_entry",&se);
    pairs.Branch("r_time",&rt);pairs.Branch("s_time",&st);pairs.Branch("r_previous_ns",&rprev);pairs.Branch("s_previous_ns",&sprev);
    pairs.Branch("r_e0_keV",&er);pairs.Branch("s_e0_keV",&es);pairs.Branch("r_old_keV",&oldr);pairs.Branch("s_old_keV",&olds);
    pairs.Branch("dt_ns",&dt);pairs.Branch("r_minor_fraction",&mr);pairs.Branch("s_minor_fraction",&ms);
    pairs.Branch("r_adjacent",&adjacent_r);pairs.Branch("s_adjacent",&adjacent_s);
    json segments=json::array(),channel_records=json::array();std::map<int,int> indices;
    for(const auto& kv:pars.channels) {
      const auto& c=kv.second;if(c.kind<3)continue;
      const int index=c.kind==3?c.ring-1:24+c.sector-1;indices[kv.first]=index;
      channel_records.push_back({{"index",index},{"key",kv.first},{"name",(c.kind==3?"ring":"sector")+rn(c.kind==3?c.ring:c.sector).substr(2)},
        {"calibration",c.cal},{"hardware_gain",c.gain},{"old_att",c.att},{"old_p0",c.fb0},{"old_p1",c.fb1}});
    }
    std::array<Long64_t,56> recorded{},selected{};Long64_t processed=0,accepted=0,pair_count=0;
    Digest sampled; // Hash exactly the decoded branch values that were inspected.
    for(block=0;block<o.blocks;++block) {
      const Long64_t low=total*block/o.blocks,high=total*(block+1)/o.blocks;
      const Long64_t count=wanted*(block+1)/o.blocks-wanted*block/o.blocks;
      const Long64_t first=low+(high-low-count)/2,last=first+count;
      tree->SetCacheEntryRange(first,last);std::vector<Hit> hits;std::map<int,Long64_t> previous;
      Long64_t native_min=std::numeric_limits<Long64_t>::max(),native_max=0,segment_s3=0;
      for(Long64_t entry=first;entry<last;++entry) {
        require(tree->GetEntry(entry)>0,"Failed decoded read");++processed;
        sampled.add(entry);sampled.add(cid);sampled.add(sid);sampled.add(ch);sampled.add(sr);sampled.add(adc);sampled.add(ts);sampled.add(pileup);sampled.add(outofr);
        if(cid!=1)continue;
        const int key=10000+100*sid+ch;auto found=indices.find(key);if(found==indices.end())continue;
        const int index=found->second;const auto& c=pars.channels.at(key);
        const int step=sr==250?8:sr==100?10:0;require(step && ts>=0 && ts<std::numeric_limits<Long64_t>::max()/step,"Invalid native time");
        Hit h;h.source=entry;h.native=ts*step;h.channel=index;h.adc=adc;h.flags=(pileup?1:0)|(outofr?2:0)|(adc==65535?4:0);
        const double x=adc;h.energy=(c.cal[0]+c.cal[1]*x+c.cal[2]*x*x)*c.gain*c.kg;
        h.old=c.energy(x,true);h.pre=c.energy(x,false);
        const double correction=timing.TotalParticleOffsetFromAdc(key,adc);
        require(std::isfinite(correction) && std::isfinite(h.energy) && std::abs(correction)<o.reorder_ns/2.,"Invalid calibration/time");
        h.time=Long64_t(double(h.native)+correction);
        native_min=std::min(native_min,h.native);native_max=std::max(native_max,h.native);++recorded[index];++segment_s3;
        hits.push_back(h);
      }
      std::sort(hits.begin(),hits.end(),[](const Hit& a,const Hit& b){return a.native<b.native || (a.native==b.native && a.source<b.source);});
      std::array<std::vector<Hit>,2> faces;
      for(auto& h:hits) {
        auto it=previous.find(h.channel);h.previous=it==previous.end()?-1:h.native-it->second;previous[h.channel]=h.native;
        if(h.pre<=o.min_energy_keV)continue;
        faces[h.channel<24?0:1].push_back(h);++selected[h.channel];++accepted;
      }
      for(auto& face:faces)std::sort(face.begin(),face.end(),[](const Hit& a,const Hit& b){return a.time<b.time || (a.time==b.time && a.source<b.source);});
      auto rings=group_face(faces[0],o.cluster_gap_ns,24),sectors=group_face(faces[1],o.cluster_gap_ns,32);
      Long64_t segment_pairs=0;size_t j0=0;
      // A 20 us boundary guard also protects the group definition and sidebands.
      for(const auto& rg:rings) {
        const auto& rh=rg.main;if(rh.native-native_min<o.reorder_ns || native_max-rh.native<o.reorder_ns)continue;
        while(j0<sectors.size() && sectors[j0].main.time<rh.time-1000)++j0;
        for(size_t j=j0;j<sectors.size() && sectors[j].main.time<=rh.time+1000;++j) {
          const auto& sg=sectors[j];const auto& sh=sg.main;
          if(sh.native-native_min<o.reorder_ns || native_max-sh.native<o.reorder_ns)continue;
          r=rh.channel;s=sh.channel;nr=rg.count;ns=sg.count;flags=rh.flags|(sh.flags<<4);group_flags=rg.flags|(sg.flags<<4);
          ra=rh.adc;sa=sh.adc;re=rh.source;se=sh.source;rt=rh.time;st=sh.time;rprev=rh.previous;sprev=sh.previous;
          er=rh.energy;es=sh.energy;oldr=rh.old;olds=sh.old;dt=double(st-rt);mr=rg.minor;ms=sg.minor;adjacent_r=rg.adjacent;adjacent_s=sg.adjacent;
          require(pairs.Fill()>0,"Pair write failed");++segment_pairs;++pair_count;
        }
      }
      segments.push_back({{"block",block},{"read_interval",{first,last}},{"source_interval",{low,high}},{"entries",count},{"s3_hits",segment_s3},{"pairs",segment_pairs},
        {"native_min_ns",native_min},{"native_max_ns",native_max},{"recorded_s3_rate_hz",native_max>native_min?double(segment_s3)*1e9/(native_max-native_min):0.}});
      std::cout<<"progress run="<<rn(o.run)<<" block="<<block+1<<" entries="<<processed<<" pairs="<<pair_count<<" seconds="<<std::chrono::duration<double>(Clock::now()-start).count()<<std::endl;
    }
    require(processed==wanted,"Incomplete sample");pars.verify();require(before.same(Stat(input)),"Input changed");
    json hashes=json::object();for(auto it=pars.snapshots.begin();it!=pars.snapshots.end();++it)hashes[it.key()]=it.value().at("sha256");
    json metadata={{"schema","sm150_s3_joint_sample_v1"},{"status","complete"},{"run",rn(o.run)},{"period",pars.group},
      {"input",input.string()},{"input_bytes",before.value.st_size},{"input_mtime_seconds",before.value.st_mtim.tv_sec},{"input_mtime_nanoseconds",before.value.st_mtim.tv_nsec},{"input_uuid",source->GetUUID().AsString()},
      {"input_entries",total},{"sample_entries",processed},{"sample_fraction",double(processed)/total},{"full_input_read",processed==total},
      {"sampled_columns_sha256",sampled.finish()},{"channels",channel_records},{"parameter_hashes",hashes},{"timing_source_sha256",timing.source_sha256()},
      {"timing_profile",timing.profile()},{"segments",segments},{"recorded_s3_hits",recorded},{"selected_s3_hits",selected},{"pairs",pair_count},
      {"min_old_pre_fb_energy_keV",o.min_energy_keV},{"face_cluster_gap_ns",o.cluster_gap_ns},{"pair_half_range_ns",1000},{"boundary_guard_ns",o.reorder_ns},
      {"energy_definition","E0 excludes attenuation and front/back; old energy is an explicit baseline, not input to E0"},
      {"selection","Per-face temporal groups; keep largest old pre-FB energy hit; no addback or front/back energy gate. Multiplicity/minor energy retained."},
      {"previous_hit_definition","Previous recorded same-channel hit in native-time order, including quality/energy rejects; -1 means censored at segment start"},
      {"elapsed_seconds",std::chrono::duration<double>(Clock::now()-start).count()},
      {"program_sha256",hash(read_text(fs::canonical("/proc/self/exe")))}};
    target->cd();require(pairs.Write()>0 && TObjString(metadata.dump().c_str()).Write("metadata")>0,"Output write failed");pairs.SetDirectory(nullptr);target->Close();
    std::unique_ptr<TFile> check(TFile::Open(partial.c_str(),"READ"));require(check && !check->IsZombie(),"Cannot reopen sample");
    auto* pt=dynamic_cast<TTree*>(check->Get("pairs"));require(pt && pt->GetEntries()==pair_count,"Pair readback failed");check->Close();
    require(::link(partial.c_str(),output.c_str())==0 && ::unlink(partial.c_str())==0,"Cannot publish sample");
    std::ofstream receipt(output.string()+".json");receipt<<metadata.dump(2)<<'\n';receipt.close();require(bool(receipt),"Cannot write receipt");
    std::cout<<"complete run="<<rn(o.run)<<" pairs="<<pair_count<<" output="<<output<<std::endl;return 0;
  }catch(const std::exception& e){std::cerr<<"build_s3 sample: "<<e.what()<<'\n';return 1;}
}
