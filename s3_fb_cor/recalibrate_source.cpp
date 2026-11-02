// Recompute S3 source energies from saved ADC values without changing pair identity.
// The exact re-expression of the adopted A/B/q is a control, not a new fit.
#include "sm150_event_calibration.h"
#include <TFile.h>
#include <TObjString.h>
#include <TTree.h>
#include <array>
#include <chrono>
#include <iostream>
#include <memory>
#include <unistd.h>
using namespace sm150::diagnostic;
namespace {
std::string file_hash(const fs::path& path) {
  std::ifstream in(path,std::ios::binary);require(bool(in),"Cannot hash input");
  Digest digest;std::array<char,1<<20> bytes{};
  while(in){in.read(bytes.data(),bytes.size());digest.add(bytes.data(),in.gcount());}
  require(in.eof(),"Input hash read failure");return digest.finish();
}
}
int main(int argc,char** argv) {
 try {
  fs::path input,output,source_path,baseline_path;std::string expected;
  for(int i=1;i<argc;++i) {
    const std::string arg=argv[i];auto value=[&](){require(++i<argc,"Missing value");return std::string(argv[i]);};
    if(arg=="--input")input=value();else if(arg=="--output")output=value();
    else if(arg=="--source-calibration")source_path=value();else if(arg=="--baseline")baseline_path=value();
    else if(arg=="--expected-sha256")expected=value();else throw std::runtime_error("Unknown option: "+arg);
  }
  require(!input.empty()&&!output.empty()&&!source_path.empty()&&!baseline_path.empty()&&expected.size()==64,"Missing required arguments");
  require(!fs::exists(output)&&!fs::exists(output.string()+".partial"),"Output exists");
  const auto start=std::chrono::steady_clock::now();const Stat before(input);
  const std::string input_hash=file_hash(input);require(input_hash==expected,"Frozen pair-input SHA mismatch");
  std::unique_ptr<TFile> src(TFile::Open(input.c_str(),"READ"));require(src&&!src->IsZombie(),"Cannot open pair input");
  auto* t=dynamic_cast<TTree*>(src->Get("pairs"));auto* obj=dynamic_cast<TObjString*>(src->Get("metadata"));require(t&&obj,"Missing pair tree/metadata");
  const json original=json::parse(obj->GetString().Data());json meta=original;
  require(original.at("schema")=="sm150_s3_joint_sample_v1","Wrong input schema");
  const std::string run=original.at("run"),period=original.at("period");const int active=period=="low"?24:8;
  const json source=json::parse(read_text(source_path)),formal=json::parse(read_text(baseline_path)).at("runs").at(run);
  require(formal.at("period")==period&&formal.at("active_rings")==active,"Baseline scope mismatch");
  std::map<int,json> source_by_key;
  for(const auto& c:source.at("channels"))if(c.at("detector")=="S3 ring"||c.at("detector")=="S3 sector") {
    const int key=period=="low"?c.at("low_beam_key").get<int>():c.at("key").get<int>();
    require(c.at("window_checked").get<bool>()&&c.at("candidate_available").get<bool>(),"Unsupported new source coefficient");
    require(source_by_key.emplace(key,c).second,"Duplicate source hardware key");
  }
  require(source_by_key.size()==56,"Incomplete new S3 source coefficients");
  std::array<double,56> c0{},c1{},g{},old0{},old1{},a{},b{},olda{},oldb{};
  std::array<bool,56> seen{};json transforms=json::array();
  for(auto& c:meta["channels"]) {
    const int i=c.at("index"),key=c.at("key");require(i>=0&&i<56&&!seen[i],"Invalid channel index");seen[i]=true;
    const auto& n=source_by_key.at(key);old0[i]=c.at("calibration")[0];old1[i]=c.at("calibration")[1];g[i]=c.at("hardware_gain");
    require(double(c.at("calibration")[2])==0.&&double(n.at("c2"))==0.,"Nonlinear source calibration needs a different rebase");
    require(std::abs(old0[i]-double(n.at("old_c0")))<1e-7&&std::abs(old1[i]-double(n.at("old_c1")))<1e-9,"Source calibration provenance mismatch");
    c0[i]=n.at("c0");c1[i]=n.at("c1");require(c1[i]>0&&old1[i]>0&&g[i]>0,"Invalid gain");
    olda[i]=formal.at("A_keV")[i];oldb[i]=formal.at("B")[i];
    const double slope=c1[i]/old1[i],offset=g[i]*(c0[i]-slope*old0[i]);
    b[i]=oldb[i]/slope;a[i]=olda[i]-b[i]*offset;
    transforms.push_back({{"index",i},{"key",key},{"channel",c.at("name")},{"production",i<active||i>=24},
      {"source_coefficients",{c0[i],c1[i],0.}},{"old_source_coefficients",{old0[i],old1[i],0.}},
      {"hardware_gain",g[i]},{"E0_new_offset_keV",offset},{"E0_new_slope",slope},
      {"A_keV",a[i]},{"B",b[i]},{"old_A_keV",olda[i]},{"old_B",oldb[i]}});
    c["calibration"]={c0[i],c1[i],0.};c["old_p0"]=a[i];c["old_p1"]=1.;c["old_att"]=b[i];
  }
  require(std::all_of(seen.begin(),seen.end(),[](bool v){return v;}),"Incomplete metadata channel map");
  Int_t block=0,r=0,s=0,nr=0,ns=0,flags=0,r_adc=0,s_adc=0;
  Long64_t r_entry=0,s_entry=0;Double_t dt=0,old_re=0,old_se=0,re=0,se=0,rb=0,sb=0;
  t->SetBranchStatus("*",false);
  auto address=[&](const char* name,auto* ptr){require(t->GetBranch(name),std::string("Missing branch ")+name);t->SetBranchStatus(name,true);require(t->SetBranchAddress(name,ptr)>=0,"Branch type mismatch");};
  address("block",&block);address("r",&r);address("s",&s);address("nr",&nr);address("ns",&ns);address("flags",&flags);
  address("r_adc",&r_adc);address("s_adc",&s_adc);address("r_entry",&r_entry);address("s_entry",&s_entry);
  address("dt_ns",&dt);address("r_e0_keV",&old_re);address("s_e0_keV",&old_se);
  fs::create_directories(fs::absolute(output).parent_path());const fs::path partial=output.string()+".partial";
  std::unique_ptr<TFile> dest(TFile::Open(partial.c_str(),"CREATE","S3 source recalibration from ADC",404));require(dest&&!dest->IsZombie(),"Cannot create output");
  TTree out("pairs","Clean production S3 pairs with new source calibration");
  out.Branch("block",&block);out.Branch("r",&r);out.Branch("s",&s);out.Branch("nr",&nr);out.Branch("ns",&ns);out.Branch("flags",&flags);
  out.Branch("r_adc",&r_adc);out.Branch("s_adc",&s_adc);out.Branch("r_entry",&r_entry);out.Branch("s_entry",&s_entry);
  out.Branch("dt_ns",&dt);out.Branch("r_e0_keV",&re);out.Branch("s_e0_keV",&se);out.Branch("r_old_keV",&rb);out.Branch("s_old_keV",&sb);
  uint64_t input_count=t->GetEntries(),clean=0,outer=0,prompt=0,random=0;double old_error=0,equiv_error=0;
  for(Long64_t entry=0;entry<t->GetEntries();++entry) {
    require(t->GetEntry(entry)>0,"Failed pair read");require(r>=0&&r<24&&s>=24&&s<56,"Invalid channel labels");
    const double er=g[r]*(old0[r]+old1[r]*r_adc),es=g[s]*(old0[s]+old1[s]*s_adc);
    old_error=std::max({old_error,std::abs(er-old_re),std::abs(es-old_se)});
    require(std::isfinite(er)&&std::isfinite(es)&&std::isfinite(dt)&&old_error<1e-6,"Original ADC energy mismatch");
    if(flags||nr!=1||ns!=1)continue;
    ++clean;
    if(r>=active){++outer;continue;}
    const bool p=dt>=-80&&dt<80,rand=(dt>=-880&&dt< -720)||(dt>=720&&dt<880);if(!p&&!rand)continue;
    re=g[r]*(c0[r]+c1[r]*r_adc);se=g[s]*(c0[s]+c1[s]*s_adc);
    rb=olda[r]+oldb[r]*er;sb=olda[s]+oldb[s]*es;
    equiv_error=std::max({equiv_error,std::abs(a[r]+b[r]*re-rb),std::abs(a[s]+b[s]*se-sb)});
    require(equiv_error<1e-6,"Equivalent transformation mismatch");
    require(out.Fill()>0,"Cannot fill output");if(p)++prompt;else ++random;
  }
  require(before.same(Stat(input)),"Input changed during read");require(uint64_t(out.GetEntries())==prompt+random,"Output count mismatch");
  meta["source_recalibration"]={{"schema","sm150_s3_source_reexpression_v1"},{"input",fs::absolute(input).string()},
    {"input_sha256",input_hash},{"source_calibration",fs::absolute(source_path).string()},{"source_sha256",hash(read_text(source_path))},
    {"baseline",fs::absolute(baseline_path).string()},{"baseline_sha256",hash(read_text(baseline_path))},{"original_metadata",original},
    {"channels",transforms},{"baseline_q",formal.at("q_position")},{"input_pairs",input_count},{"clean_pairs",clean},{"excluded_outer_clean_pairs",outer},
    {"prompt",prompt},{"random",random},{"old_ADC_energy_max_error_keV",old_error},{"equivalent_energy_max_error_keV",equiv_error},
    {"branch_definition","r/s_e0_keV=new ADC source energy; r/s_old_keV=current formal A+B*old_E0 before q, used only as frozen comparison and fit coordinate"},
    {"selection","flags=0,nr=ns=1, production rings, original prompt/random times; no additional energy or FB gate"},
    {"program_sha256",file_hash(fs::canonical("/proc/self/exe"))},
    {"elapsed_seconds",std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()}};
  out.Write();TObjString metadata(meta.dump().c_str());metadata.Write("metadata");out.SetDirectory(nullptr);dest->Close();src->Close();
  require(::link(partial.c_str(),output.c_str())==0&&::unlink(partial.c_str())==0,"Cannot publish without overwriting");
  std::ofstream receipt(output.string()+".json");receipt<<meta.dump(2)<<'\n';receipt.close();require(bool(receipt),"Cannot write receipt");
  std::cout<<run<<" prepared "<<prompt+random<<" pairs; max energy error="<<equiv_error<<" keV\n";
  return 0;
 }catch(const std::exception& e){std::cerr<<"S3 source recalibration: "<<e.what()<<'\n';return 1;}
}
