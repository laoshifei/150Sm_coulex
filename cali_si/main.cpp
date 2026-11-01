#include <TROOT.h>
#include "get_par.cpp"
#include "compute_average.cpp"
int main(int argc,char** argv) {
  try {
    gROOT->SetBatch(true);sm150::sourcecal::Options options(argc,argv);
    if(options.help){sm150::sourcecal::help("cali_si");std::cout<<"Additional commands: average --input calibration.json --keys K1,K2 --output NEW_DIR; remap-low --input RUN0457.json --output NEW_DIR\n";return 0;}
    if(options.command=="source")return get_par_all(options);
    if(options.command=="average")return compute_average(options);
    if(options.command=="remap-low")return remap_low(options);
    throw std::runtime_error("Unknown command: "+options.command);
  } catch(const std::exception& e){std::cerr<<"cali_si: "<<e.what()<<'\n';return 1;}
}
