#include "build.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
int main(int argc,char** argv) {
  try {
    BuildOptions o;
    for(int i=1;i<argc;++i) {
      std::string a=argv[i];
      if(a=="--help" || a=="-h") {
        std::cout<<"./build RUN [--input decoded.root] [--output s3.root] [--project DIR]\n"
          "  --profile-dir DIR  Complete timing profile; default active.txt\n"
          "  --entries N        Prefix limit for tests; omitted means full input\n"
          "  --sample-entries N Distributed joint-calibration pairs; 0 means all\n"
          "  --gap-ns N         Temporal connected-component gap (default 100)\n"
          "  --reorder-ns N     Bounded reorder buffer (default 20000)\n"
          "  --min-energy-kev E Pre-normalization threshold (default 2500)\n"
          "  --blocks N         Contiguous fit/check blocks (default 40)\n"
          "No addback or energy-agreement gate; existing outputs are protected.\n";
        return 0;
      }
      auto value=[&](){if(++i>=argc)throw std::invalid_argument("Missing value for "+a);return std::string(argv[i]);};
      auto integer=[](const std::string& s){size_t n=0;auto v=std::stoll(s,&n);if(n!=s.size())throw std::invalid_argument("Invalid integer");return v;};
      if(a=="--input")o.input=value();else if(a=="--output")o.output=value();
      else if(a=="--project")o.project=value();else if(a=="--profile-dir")o.profile_dir=value();
      else if(a=="--entries")o.entries=integer(value());else if(a=="--gap-ns")o.cluster_gap_ns=integer(value());
      else if(a=="--sample-entries")o.sample_entries=integer(value());
      else if(a=="--reorder-ns")o.reorder_ns=integer(value());else if(a=="--blocks")o.blocks=integer(value());
      else if(a=="--min-energy-kev") {auto s=value();size_t n=0;o.min_energy_keV=std::stod(s,&n);if(n!=s.size())throw std::invalid_argument("Invalid energy");}
      else if(!a.empty() && a[0]!='-' && !o.run)o.run=integer(a);
      else throw std::invalid_argument("Unknown argument: "+a);
    }
    if(o.run<1 || o.run>9999 || o.entries<0 || o.cluster_gap_ns<1 || o.cluster_gap_ns>2000 ||
       o.reorder_ns<2*o.cluster_gap_ns || o.blocks<2 || o.blocks>1000 || !std::isfinite(o.min_energy_keV) || o.min_energy_keV<0.)
      throw std::invalid_argument("Invalid settings; use --help");
    if(o.sample_entries < -1 || (o.sample_entries>=0 && o.entries))
      throw std::invalid_argument("Distributed sampling and prefix limits cannot be combined");
    return o.sample_entries>=0?build(o).Sample():build(o).Process();
  }catch(const std::exception& e){std::cerr<<"build_s3: "<<e.what()<<'\n';return 1;}
}
