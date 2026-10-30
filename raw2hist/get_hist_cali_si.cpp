#include "run_histogram.h"
// This compatibility entry delegates to raw2hist/histogram.
void get_hist_cali_si(TString filename, int run) {
  raw2hist_macro::run(run, filename.Data(), "both", "particle", "decoded");
}
