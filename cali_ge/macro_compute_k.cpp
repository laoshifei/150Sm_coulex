#include "../include/sm150_source_calibration_cli.h"
int macro_compute_k(int run=0,const char* input="",const char* seed="",const char* output="") {
  if(run<=0||!*input||!*seed||!*output) {
    std::cerr<<"Usage: macro_compute_k(run, modern_base_histogram, base_seed_table, new_output_directory)\n";return 1;
  }
  return sm150_source_cli::run(__FILE__,{"gain","--run",std::to_string(run),"--input",input,"--seed",seed,"--output",output});
}
