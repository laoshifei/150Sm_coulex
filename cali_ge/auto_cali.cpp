#ifndef SM150_GE_AUTO_CALI_CPP
#define SM150_GE_AUTO_CALI_CPP
#include "../include/sm150_source_calibration.h"
namespace sm150 { namespace sourcecal {
// The line policy is explicit and independent of ROOT session globals.
inline json ge_channel_config(int, const json& config) { return config; }
}}
#endif
