#include "auto_cali.cpp"
int get_par_all(const sm150::sourcecal::Options& options) {
  return sm150::sourcecal::run_source(options,0,sm150::sourcecal::ge_channel_config);
}
