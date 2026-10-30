#include "run_histogram.h"
// Historical floating-energy input; no second calibration is applied.
void get_hist(int run) {
  auto input=raw2hist_macro::project()/"rootfile"/TString::Format("data%04d.root",run).Data();
  raw2hist_macro::run(run,input.string(),"calibrated","all","calibrated");
}
