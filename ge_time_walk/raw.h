#ifndef SM150_OFFSET_RAW_H
#define SM150_OFFSET_RAW_H
#include <TTree.h>
// evte 始终是原始 ADC。刻度能量另用 double，不写回 16 位变量。
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
