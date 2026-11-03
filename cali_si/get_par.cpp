#include "auto_cali.cpp"
#include "joint_source.h"
int get_par_all(const sm150::sourcecal::Options& options) {
  return sm150::sourcecal::run_source(options,1,sm150::sourcecal::si_channel_config,sm150::sourcecal::joint_alpha_channel);
}
