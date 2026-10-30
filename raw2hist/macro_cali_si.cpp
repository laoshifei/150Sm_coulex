#include "run_histogram.h"
// The default input is decoded C2 on E; an explicit path can also select C1.
void macro_cali_si(int run, const char* input="") {
  raw2hist_macro::run(run,input,"both","particle","decoded");
}
