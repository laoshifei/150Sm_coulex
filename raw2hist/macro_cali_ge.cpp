#include "run_histogram.h"
// The default input is decoded C2 on E; an explicit path can also select C1.
void macro_cali_ge(int run, const char* input="") {
  raw2hist_macro::run(run,input,"both","ge","decoded");
}
