#include <TROOT.h>
#include "get_par.cpp"
#include "compute_k.cpp"
int main(int argc,char** argv) {
  try {
    gROOT->SetBatch(true);sm150::sourcecal::Options options(argc,argv);
    if(options.help){sm150::sourcecal::help("cali_ge");std::cout<<"Use 'gain' with modern energy_base histograms for 511-keV run gains.\n";return 0;}
    if(options.command=="source")return get_par_all(options);
    if(options.command=="gain")return compute_k_all(options);
    throw std::runtime_error("Unknown command: "+options.command);
  } catch(const std::exception& e){std::cerr<<"cali_ge: "<<e.what()<<'\n';return 1;}
}
