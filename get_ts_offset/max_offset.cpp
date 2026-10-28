#include "offset_fit.h"
// 最高 bin 用 GetBinCenter；没有伪造的拟合 chi2 或统计误差。
int max_offset(TString filename,int run,TString outdir="offset_peak_candidates",
               TString selection="isolated",int min_entries=100,double halfwidth_ns=100.,
               TString scope="legacy") {
  return sm150_offset::export_offsets(filename,run,outdir,false,selection.Data(),min_entries,halfwidth_ns,scope.Data());
}
