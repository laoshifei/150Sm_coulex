#include "offset_fit.h"
// Use GetBinCenter for the highest bin; do not invent a fit chi2 or statistical uncertainty.
int max_offset(TString filename,int run,TString outdir="offset_peak_candidates",
               TString selection="isolated",int min_entries=100,double halfwidth_ns=100.,
               TString scope="legacy") {
  return sm150_offset::export_offsets(filename,run,outdir,false,selection.Data(),min_entries,halfwidth_ns,scope.Data());
}
