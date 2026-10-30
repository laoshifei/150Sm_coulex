#include "s3_tree.h"
#include "sm150_event_calibration.h"
#include <set>
using sm150::diagnostic::require;
namespace {
template<class T> bool ready(const std::unique_ptr<T>& value) {
  // TTreeReader loads branches lazily; IsValid is false before the first read.
  return value->GetSetupStatus()>=0 &&
    value->ProxyRead()==ROOT::Internal::TTreeReaderValueBase::kReadSuccess;
}
}

s3_tree::s3_tree(TTree* tree,const std::string& format):reader_(tree),modern_(format=="s3_events") {
  require(tree,"Missing input tree");
  tree->SetBranchStatus("*",false);
  tree->SetCacheSize(8*1024*1024);
  auto enable=[&](const char* name) {
    require(tree->GetBranch(name),std::string("Missing branch: ")+name);
    tree->SetBranchStatus(name,true);tree->AddBranchToCache(name,true);
  };
  if(modern_) {
    for(auto name:{"cluster","boundary","n_ring_raw","n_sector_raw","face","id","channel","adc","flags","source_entry","time_ns","energy_previous_keV"})enable(name);
    nr_=std::make_unique<TTreeReaderValue<Int_t>>(reader_,"n_ring_raw");
    ns_=std::make_unique<TTreeReaderValue<Int_t>>(reader_,"n_sector_raw");
    cluster_=std::make_unique<TTreeReaderValue<Long64_t>>(reader_,"cluster");
    boundary_=std::make_unique<TTreeReaderValue<Bool_t>>(reader_,"boundary");
    face_=std::make_unique<TTreeReaderValue<std::vector<int>>>(reader_,"face");
    id_=std::make_unique<TTreeReaderValue<std::vector<int>>>(reader_,"id");
    key_=std::make_unique<TTreeReaderValue<std::vector<int>>>(reader_,"channel");
    adc_=std::make_unique<TTreeReaderValue<std::vector<int>>>(reader_,"adc");
    flags_=std::make_unique<TTreeReaderValue<std::vector<int>>>(reader_,"flags");
    source_=std::make_unique<TTreeReaderValue<std::vector<Long64_t>>>(reader_,"source_entry");
    time_=std::make_unique<TTreeReaderValue<std::vector<Long64_t>>>(reader_,"time_ns");
    energy_=std::make_unique<TTreeReaderValue<std::vector<double>>>(reader_,"energy_previous_keV");
  } else {
    for(auto name:{"n_s3_ring","n_s3_sector","s3_ring_id","s3_sector_id","s3_ring_energy","s3_sector_energy","s3_ring_ts","s3_sector_ts"})enable(name);
    nr_=std::make_unique<TTreeReaderValue<Int_t>>(reader_,"n_s3_ring");
    ns_=std::make_unique<TTreeReaderValue<Int_t>>(reader_,"n_s3_sector");
    ring_id_=std::make_unique<TTreeReaderArray<Short_t>>(reader_,"s3_ring_id");
    sector_id_=std::make_unique<TTreeReaderArray<Short_t>>(reader_,"s3_sector_id");
    ring_energy_=std::make_unique<TTreeReaderArray<Double_t>>(reader_,"s3_ring_energy");
    sector_energy_=std::make_unique<TTreeReaderArray<Double_t>>(reader_,"s3_sector_energy");
    ring_time_=std::make_unique<TTreeReaderArray<Long64_t>>(reader_,"s3_ring_ts");
    sector_time_=std::make_unique<TTreeReaderArray<Long64_t>>(reader_,"s3_sector_ts");
  }
  tree->StopCacheLearningPhase();
}

S3Event s3_tree::Read(Long64_t entry) {
  const bool ok=entry==previous_+1?reader_.Next():reader_.SetEntry(entry)==TTreeReader::kEntryValid;
  require(ok,"Cannot read source entry "+std::to_string(entry));previous_=entry;
  require(ready(nr_)&&ready(ns_),"Invalid multiplicity branch types");
  const int nr=**nr_,ns=**ns_;require(nr>=0&&ns>=0,"Negative multiplicity");
  S3Event out;out.entry=entry;out.cluster=entry;
  if(modern_) {
    require(ready(cluster_)&&ready(boundary_)&&ready(face_)&&ready(id_)&&ready(key_)&&ready(adc_)&&ready(flags_)&&ready(source_)&&ready(time_)&&ready(energy_),"Invalid s3_events branch types");
    const auto n=(**face_).size();
    require(n==size_t(nr+ns)&&(**id_).size()==n&&(**key_).size()==n&&(**adc_).size()==n&&(**flags_).size()==n&&(**source_).size()==n&&(**time_).size()==n&&(**energy_).size()==n,"Inconsistent vector sizes");
    out.cluster=**cluster_;out.boundary=**boundary_;std::set<Long64_t> seen;
    for(size_t i=0;i<n;++i) {
      int face=(**face_)[i];require(face==0||face==1,"Invalid S3 face");
      S3Hit hit;hit.id=(**id_)[i];hit.key=(**key_)[i];hit.adc=(**adc_)[i];hit.flags=(**flags_)[i];
      hit.source=(**source_)[i];hit.time=(**time_)[i];hit.energy=(**energy_)[i];
      require(hit.source>=0&&seen.insert(hit.source).second,"Duplicate source hit in cluster");
      (face?out.sectors:out.rings).push_back(hit);
    }
  } else {
    // Array implementations perform their own read when GetSize/At is called.
    require(ring_id_->GetSetupStatus()>=0&&sector_id_->GetSetupStatus()>=0&&ring_energy_->GetSetupStatus()>=0&&sector_energy_->GetSetupStatus()>=0&&ring_time_->GetSetupStatus()>=0&&sector_time_->GetSetupStatus()>=0,"Invalid legacy branch types");
    require(ring_id_->GetSize()==size_t(nr)&&ring_energy_->GetSize()==size_t(nr)&&ring_time_->GetSize()==size_t(nr)&&sector_id_->GetSize()==size_t(ns)&&sector_energy_->GetSize()==size_t(ns)&&sector_time_->GetSize()==size_t(ns),"Inconsistent legacy arrays");
    for(int face=0;face<2;++face)for(int i=0;i<(face?ns:nr);++i) {
      S3Hit hit;hit.id=face?(*sector_id_)[i]:(*ring_id_)[i];
      hit.energy=face?(*sector_energy_)[i]:(*ring_energy_)[i];hit.time=face?(*sector_time_)[i]:(*ring_time_)[i];
      (face?out.sectors:out.rings).push_back(hit);
    }
  }
  require(out.rings.size()==size_t(nr)&&out.sectors.size()==size_t(ns),"Face count mismatch");
  for(int face=0;face<2;++face)for(const auto& hit:face?out.sectors:out.rings)
    require(hit.id>=1&&hit.id<=(face?32:24)&&std::isfinite(hit.energy)&&std::isfinite(hit.time),"Invalid S3 hit");
  return out;
}
