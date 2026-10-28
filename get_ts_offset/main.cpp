#include "analysis.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace {
uint64_t integer(const std::string& value) {
  if (value.empty() || value.find_first_not_of("0123456789")!=std::string::npos)
    throw std::runtime_error("Expected non-negative integer: "+value);
  return std::stoull(value);
}
double number(const std::string& value) {
  size_t n=0; double x=std::stod(value,&n);
  if(n!=value.size() || !std::isfinite(x)) throw std::runtime_error("Invalid number: "+value);
  return x;
}
void help() {
  std::cout << "Usage: analysis RUN --input C2.root --output candidate.root [options]\n"
    "  --project DIR              project containing pars (auto-detected)\n"
    "  --mode raw|residual        native timestamps / current complete correction once\n"
    "  --profile NAME            explicit residual profile; default active.txt\n"
    "  --energy-stage final|base  final calibration / omit S3 front-back affine term\n"
    "  --ge-min/--ge-max KEV      default [40,2000)\n"
    "  --si-min/--si-max KEV      default (15000,100000)\n"
    "  --s3-min/--s3-max KEV      default (25000,100000)\n"
    "  --window NS --bin NS       strict |dt|<window; default 1000 / 2\n"
    "  --reference KEY           Ge reference; default 212, normal CFD only\n"
    "  --begin N --entries N      contiguous interval; entries=0 to EOF\n"
    "  --chunk N                 default 100000, bounded rolling buffer\n"
    "  --max-buffer-hits N        fail explicitly if buffer exceeds budget\n"
    "  --write-pairs             source entries, dt, energies and multiplicities\n"
    "  --particle-calibration    compact per-particle Ge / S3-face diagnostics\n"
    "  --particle-network        all Si-Si and Si-S3 pairs, no Ge requirement\n"
    "  --particle-full           full-statistics particle histograms, no pair tree\n"
    "  --inspect-parameters      print calibration; no ROOT event reading\n"
    "Existing outputs are not overwritten. No parameters are installed.\n";
}
}
int main(int argc, char** argv) {
  try {
    if(argc<2 || std::string(argv[1])=="--help") { help(); return argc<2?1:0; }
    OffsetOptions o;
    auto run=integer(argv[1]);
    if(run==0 || run>9999) throw std::runtime_error("Invalid run");
    o.run=int(run);
    for(int i=2;i<argc;++i) {
      const std::string k=argv[i];
      if(k=="--write-pairs") { o.write_pairs=true; continue; }
      if(k=="--inspect-parameters") { o.inspect_parameters=true; continue; }
      if(k=="--particle-calibration") { o.particle_calibration=true; continue; }
      if(k=="--particle-network") { o.particle_network=true; continue; }
      if(k=="--particle-full") { o.particle_full=true; continue; }
      if(k=="--help") { help(); return 0; }
      if(i+1>=argc) throw std::runtime_error("Missing value for "+k);
      const std::string v=argv[++i];
      if(k=="--project") o.project=v;
      else if(k=="--input") o.input=v;
      else if(k=="--output") o.output=v;
      else if(k=="--mode") o.mode=v;
      else if(k=="--profile") o.profile=v;
      else if(k=="--energy-stage") o.energy_stage=v;
      else if(k=="--reference") { auto n=integer(v); if(n>9999) throw std::runtime_error("Invalid Ge key"); o.ge_reference=int(n); }
      else if(k=="--window") o.window_ns=number(v);
      else if(k=="--bin") o.bin_ns=number(v);
      else if(k=="--ge-min") o.ge_min=number(v);
      else if(k=="--ge-max") o.ge_max=number(v);
      else if(k=="--si-min") o.si_min=number(v);
      else if(k=="--si-max") o.si_max=number(v);
      else if(k=="--s3-min") o.s3_min=number(v);
      else if(k=="--s3-max") o.s3_max=number(v);
      else if(k=="--begin") o.begin=integer(v);
      else if(k=="--entries") o.entries=integer(v);
      else if(k=="--chunk") o.chunk=integer(v);
      else if(k=="--max-buffer-hits") o.max_buffer_hits=integer(v);
      else throw std::runtime_error("Unknown option: "+k);
    }
    if(o.project.empty()) {
      for(const auto& p : {std::filesystem::current_path(),std::filesystem::current_path().parent_path(),
           std::filesystem::canonical("/proc/self/exe").parent_path().parent_path()})
        if(std::filesystem::is_directory(p/"pars")) { o.project=p; break; }
    }
    if(o.project.empty()) throw std::runtime_error("Use --project to locate pars");
    analysis(o).Process();
    return 0;
  } catch(const std::exception& e) { std::cerr<<"ERROR: "<<e.what()<<'\n'; return 1; }
}
