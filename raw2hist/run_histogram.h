#ifndef SM150_RAW2HIST_RUN_H
#define SM150_RAW2HIST_RUN_H
// Compatibility macros delegate to the single compiled implementation.
#include <TString.h>
#include <TSystem.h>
#include <filesystem>
#include <stdexcept>
#include <string>
namespace raw2hist_macro {
inline std::filesystem::path project() {
  auto p=std::filesystem::current_path();
  while(!p.empty()) {
    if(std::filesystem::is_directory(p/"raw2hist") && std::filesystem::is_directory(p/"pars"))return p;
    if(p==p.parent_path())break;
    p=p.parent_path();
  }
  throw std::runtime_error("Run this macro from the project or raw2hist directory");
}
inline std::string quote(const std::string& s) {
  std::string out="'";for(char c:s)out+=c=='\''?"'\\''":std::string(1,c);return out+"'";
}
inline void run(int number,const std::string& input,const std::string& mode,
                const std::string& detector="all",const std::string& format="decoded") {
  auto p=project();auto exe=p/"raw2hist/histogram";
  if(!std::filesystem::is_regular_file(exe))throw std::runtime_error("Compile first: make -C raw2hist");
  const auto tag=TString::Format("%04d",number);const std::string id=tag.Data();
  auto output=p/"rootfile/raw2hist"/("data"+id+"_"+format+"_"+detector+"_"+mode+".root");
  std::string cmd=quote(exe.string())+" --project "+quote(p.string())+" --run "+id+" --output "+quote(output.string())+
    " --mode "+mode+" --detector "+detector+" --format "+format;
  if(!input.empty())cmd+=" --input "+quote(input);
  if(gSystem->Exec(cmd.c_str())!=0)throw std::runtime_error("raw2hist failed; inspect the preceding error");
}
}
#endif
