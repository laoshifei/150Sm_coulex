#include "build.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace {
uint64_t integer(const std::string& s) {
  if(s.empty() || s.find_first_not_of("0123456789")!=std::string::npos)
    throw std::runtime_error("Expected unsigned integer: "+s);
  return std::stoull(s);
}
double number(const std::string& s) {
  size_t n=0;double x=std::stod(s,&n);
  if(n!=s.size() || !std::isfinite(x)) throw std::runtime_error("Invalid number: "+s);
  return x;
}
void help() {
  std::cout<<"Usage: build RUN --input C2.root --output walk.root [options]\n"
    "  --project DIR                 project containing energy calibration pars/\n"
    "  --profile local|active|NAME   local package (default) or installed version\n"
    "  --profile-dir DIR            explicit complete package; excludes --profile\n"
    "  --begin N --entries N          diagnostic interval; 0 entries=EOF\n"
    "  --chunk N                     default100000\n"
    "  --native-window NS            strict coarse pair window, default2000\n"
    "  --time-range NS --time-bin NS  histograms default+/-1000,4\n"
    "  --energy-bin KEV              gamma bin default10\n"
    "  --ge-min/--ge-max KEV          default[40,2000)\n"
    "  --si-min/--si-max KEV          default(15000,100000)\n"
    "  --s3-min/--s3-max KEV          default(25000,100000)\n"
    "  --s3-match-ns NS              independent ring/sector match, default80\n"
    "  --s3-match-keV KEV            front/back energy difference, default2500\n"
    "  --max-buffer-hits N           explicit memory guard\n"
    "  --write-pairs                 source entries, energies and raw CFD fields\n"
    "  --inspect-parameters          print calibration, no event reading\n"
    "All Ge-particle candidates, no nearest-time choice; normal/failed retained.\n"
    "Existing outputs refused; current/native-Ge views use identical pairs.\n";
}
}
int main(int argc,char** argv) {
  try {
    if(argc<2 || std::string(argv[1])=="--help") { help();return argc<2?1:0; }
    WalkOptions o;auto run=integer(argv[1]);
    if(!run || run>9999) throw std::runtime_error("Invalid run");
    o.run=int(run);
    for(int i=2;i<argc;++i) {
      std::string k=argv[i];
      if(k=="--write-pairs") {o.write_pairs=true;continue;}
      if(k=="--inspect-parameters") {o.inspect_parameters=true;continue;}
      if(k=="--help") {help();return 0;}
      if(i+1>=argc) throw std::runtime_error("Missing value for "+k);
      const std::string v=argv[++i];
      if(k=="--project") o.project=v;
      else if(k=="--input") o.input=v;
      else if(k=="--output") o.output=v;
      else if(k=="--profile") o.profile=v;
      else if(k=="--profile-dir") o.profile_dir=v;
      else if(k=="--begin") o.begin=integer(v);
      else if(k=="--entries") o.entries=integer(v);
      else if(k=="--chunk") o.chunk=integer(v);
      else if(k=="--max-buffer-hits") o.max_buffer_hits=integer(v);
      else if(k=="--native-window") o.native_window_ns=number(v);
      else if(k=="--time-range") o.time_range_ns=number(v);
      else if(k=="--time-bin") o.time_bin_ns=number(v);
      else if(k=="--energy-bin") o.energy_bin_keV=number(v);
      else if(k=="--ge-min") o.ge_min=number(v);
      else if(k=="--ge-max") o.ge_max=number(v);
      else if(k=="--si-min") o.si_min=number(v);
      else if(k=="--si-max") o.si_max=number(v);
      else if(k=="--s3-min") o.s3_min=number(v);
      else if(k=="--s3-max") o.s3_max=number(v);
      else if(k=="--s3-match-ns") o.s3_match_ns=number(v);
      else if(k=="--s3-match-keV") o.s3_match_keV=number(v);
      else throw std::runtime_error("Unknown option: "+k);
    }
    if(o.project.empty()) for(const auto& p:{std::filesystem::current_path(),std::filesystem::current_path().parent_path(),std::filesystem::canonical("/proc/self/exe").parent_path().parent_path()})
      if(std::filesystem::exists(p/"pars/time_alignment/active.txt")) {o.project=p;break;}
    if(o.project.empty()) throw std::runtime_error("Use --project");
    if(!o.profile.empty() && !o.profile_dir.empty()) throw std::runtime_error("Choose --profile or --profile-dir, not both");
    build_ge_time_walk(o);return 0;
  } catch(const std::exception& e) {std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
}
