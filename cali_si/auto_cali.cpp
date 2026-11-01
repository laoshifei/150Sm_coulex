#ifndef SM150_SI_AUTO_CALI_CPP
#define SM150_SI_AUTO_CALI_CPP
#include "../include/sm150_source_calibration.h"
namespace sm150 { namespace sourcecal {
inline json si_channel_config(int key,const json& config) {
  json selected=config;
  selected["lines"]=config.at(key%10000/100<=8?"three_alpha_lines":"two_alpha_lines");
  return selected;
}
}}
#endif
