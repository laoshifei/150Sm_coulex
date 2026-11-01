#ifndef SM150_SOURCE_CALIBRATION_CLI_H
#define SM150_SOURCE_CALIBRATION_CLI_H
// Small explicit ROOT wrappers for the compiled source-calibration tools.
#include <TSystem.h>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
namespace sm150_source_cli {
inline std::string quote(const std::string& text) {
  std::string result="'";
  for(char c:text)result+=c=='\''?"'\\''":std::string(1,c);
  return result+"'";
}
inline int run(const char* macro,const std::vector<std::string>& arguments) {
  const auto binary=std::filesystem::absolute(macro).parent_path()/"calibrate";
  if(!std::filesystem::is_regular_file(binary)){std::cerr<<"Compile first: make -C "<<binary.parent_path()<<'\n';return 1;}
  std::string command=quote(binary.string());
  for(const auto& value:arguments)command+=" "+quote(value);
  return gSystem->Exec(command.c_str());
}
}
#endif
