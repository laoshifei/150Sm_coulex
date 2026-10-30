#ifndef SM150_S3_COIN_TREE_H
#define SM150_S3_COIN_TREE_H
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <memory>
#include <string>
#include <vector>

struct S3Hit {
  int id=0, key=0, adc=0, flags=0;
  Long64_t source=-1;
  double energy=0., time=0.;
};
struct S3Event {
  Long64_t entry=0, cluster=0;
  bool boundary=false;
  std::vector<S3Hit> rings, sectors;
};

// Dynamic readers prevent a legacy two-element array from truncating new hits.
class s3_tree {
 public:
  s3_tree(TTree* tree, const std::string& format);
  S3Event Read(Long64_t entry);
 private:
  TTreeReader reader_;
  bool modern_;
  Long64_t previous_=-2;
  template<class T> using Value=std::unique_ptr<TTreeReaderValue<T>>;
  template<class T> using Array=std::unique_ptr<TTreeReaderArray<T>>;
  Value<Int_t> nr_, ns_;
  Value<Long64_t> cluster_;
  Value<Bool_t> boundary_;
  Value<std::vector<int>> face_, id_, key_, adc_, flags_;
  Value<std::vector<Long64_t>> source_, time_;
  Value<std::vector<double>> energy_;
  Array<Short_t> ring_id_, sector_id_;
  Array<Double_t> ring_energy_, sector_energy_;
  Array<Long64_t> ring_time_, sector_time_;
};
#endif
