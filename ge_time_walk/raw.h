#ifndef SM150_OFFSET_RAW_H
#define SM150_OFFSET_RAW_H
#include <TTree.h>
// evte always holds raw ADC values. Store calibrated energy separately as double; never write it back to the 16-bit field.
class raw {
public:
  explicit raw(TTree* tree);
  ~raw();
  void GetEntry(Long64_t entry);
  Short_t sr{}, cid{}, sid{}, ch{}, cfd{}, cfds{};
  Bool_t pileup{}, outofr{}, cfdft{};
  UShort_t evte{};
  Long64_t ts{};
private:
  TTree* tree_;
};
#endif
