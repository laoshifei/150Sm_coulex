#include "analysis.h"
#include <iostream>
#include <stdexcept>

int main(int argc,char** argv) {
  try {
    AnalysisOptions o;
    for(int i=1;i<argc;++i) {
      const std::string arg=argv[i];
      auto value=[&]() { if(++i>=argc)throw std::runtime_error("Missing value for "+arg);return std::string(argv[i]); };
      auto real=[&]() { const auto s=value();size_t n=0;double x=std::stod(s,&n);if(n!=s.size())throw std::runtime_error("Invalid number: "+s);return x; };
      auto integer=[&]() { const auto s=value();size_t n=0;long long x=std::stoll(s,&n);if(n!=s.size())throw std::runtime_error("Invalid integer: "+s);return x; };
      if(arg=="--help") {
        std::cout<<"ana RUN --input FILE --output NEW.root [options]\n"
          <<"  --project DIR                 geometry table directory root\n"
          <<"  --format auto|s3_events|legacy explicit legacy opt-in\n"
          <<"  --sample-entries N --segments 20  distributed rows;0 means all\n"
          <<"  --s3-energy-file FILE         optional complete A/B/q candidate\n"
          <<"  --energy-min-mev 2.5 --energy-max-mev 100\n"
          <<"  --face-time-ns 80 --face-energy-mev 2.5 --particle-time-ns 80\n"
          <<"  --phi-low-deg 160 --phi-high-deg 200 --cm-low-deg 160 --cm-high-deg 200\n"
          <<"  --beam-mev 102 --deficit-mev 2 --projectile-a 32 --target-a 150\n"
          <<"Input must retain both particles. Current Ge-gated build_all files are rejected.\n";
        return 0;
      } else if(arg=="--input")o.input=value();
      else if(arg=="--output")o.output=value();
      else if(arg=="--project")o.project=value();
      else if(arg=="--format")o.format=value();
      else if(arg=="--s3-energy-file")o.energy_file=value();
      else if(arg=="--sample-entries")o.sample_entries=integer();
      else if(arg=="--segments") {auto n=integer();if(n<1||n>100000)throw std::runtime_error("Invalid segment count");o.segments=static_cast<int>(n);}
      else if(arg=="--energy-min-mev")o.energy_min=1000*real();
      else if(arg=="--energy-max-mev")o.energy_max=1000*real();
      else if(arg=="--face-time-ns")o.face_time=real();
      else if(arg=="--face-energy-mev")o.face_energy=1000*real();
      else if(arg=="--particle-time-ns")o.particle_time=real();
      else if(arg=="--phi-low-deg")o.phi_low=real();
      else if(arg=="--phi-high-deg")o.phi_high=real();
      else if(arg=="--cm-low-deg")o.cm_low=real();
      else if(arg=="--cm-high-deg")o.cm_high=real();
      else if(arg=="--beam-mev")o.beam=real();
      else if(arg=="--deficit-mev")o.deficit=real();
      else if(arg=="--projectile-a")o.projectile_mass=real();
      else if(arg=="--target-a")o.target_mass=real();
      else if(o.run==0&&arg.find_first_not_of("0123456789")==std::string::npos)o.run=std::stoi(arg);
      else throw std::runtime_error("Unknown argument: "+arg);
    }
    return analysis(o).Process();
  } catch(const std::exception& error) {
    std::cerr<<"s3_coin_ana: "<<error.what()<<'\n';return 1;
  }
}
