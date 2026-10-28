#include "raw.h"
#include <TBranch.h>
#include <TLeaf.h>
#include <stdexcept>
#include <string>
raw::raw(TTree* tree) : tree_(tree) {
  if (!tree_) throw std::runtime_error("Missing input tree");
  tree_->SetBranchStatus("*", 0);
  auto bind = [&](const char* name, const char* type, void* address) {
    auto* branch = tree_->GetBranch(name);
    auto* leaf = tree_->GetLeaf(name);
    if (!branch || !leaf || leaf->GetLeafCount() || leaf->GetLenStatic()!=1 ||
        std::string(leaf->GetTypeName()) != type)
      throw std::runtime_error(std::string("Missing/wrong scalar branch: ")+name);
    tree_->SetBranchStatus(name, 1);
    if (tree_->SetBranchAddress(name, address)<0)
      throw std::runtime_error(std::string("Cannot bind: ")+name);
  };
  bind("sr", "Short_t", &sr); bind("cid", "Short_t", &cid);
  bind("sid", "Short_t", &sid); bind("ch", "Short_t", &ch);
  bind("cfd", "Short_t", &cfd); bind("cfds", "Short_t", &cfds);
  bind("pileup", "Bool_t", &pileup); bind("outofr", "Bool_t", &outofr);
  bind("cfdft", "Bool_t", &cfdft); bind("evte", "UShort_t", &evte);
  bind("ts", "Long64_t", &ts);
}
raw::~raw() { tree_->ResetBranchAddresses(); } // 文件由调用方持有。
void raw::GetEntry(Long64_t entry) {
  if (tree_->GetEntry(entry)<=0) throw std::runtime_error("Failed reading entry "+std::to_string(entry));
}
