#include "doppler.h"
#include <iostream>

#include "set.h"

int main(int argc, char* argv[])
{
  if(argc==2 && std::string(argv[1])=="--help") {
    std::cout << "Current build: doppler --project DIR --parameters JSON --input BUILD.root --output NEW.root\n"
              << "Scan cached pairs: doppler --project DIR --scan-parameters JSON --scan-list TXT --output NEW.root\n"
              << "Historical input: doppler RUN WIN JUMP (legacy layout and selection only)\n";
    return 0;
  }
  if(argc>1 && std::string(argv[1]).rfind("--",0)==0) {
    try {
      sm150_doppler::Options o;
      for(int i=1;i<argc;i+=2) {
        if(i+1>=argc) throw std::runtime_error("Missing option value");
        const std::string key=argv[i],value=argv[i+1];
        if(key=="--project")o.project=value;
        else if(key=="--parameters")o.parameters=value;
        else if(key=="--input")o.input=value;
        else if(key=="--output")o.output=value;
        else if(key=="--scan-list")o.scan_list=value;
        else if(key=="--scan-parameters")o.scan_parameters=value;
        else throw std::runtime_error("Unknown option: "+key);
      }
      if(o.project.empty()||o.output.empty()) throw std::runtime_error("Project and output are required");
      if(!o.scan_list.empty()) return sm150_doppler::scan_pairs(o);
      return sm150_doppler::process_build(o);
    } catch(const std::exception &e) { std::cerr<<"Doppler error: "<<e.what()<<"\n"; return 1; }
  }
  if(argc!= 4){
    std::cout << "./doppler run win jump" << std::endl;
    return -1;
  }

  int run = atoi(argv[1]);
  int win = atoi(argv[2]);
  int jump = atoi(argv[3]);

#ifdef TW  
  TString fi_name = TString::Format("../rootfile/data%04d_tw_build_%dns_jump_%dns.root", run, win, jump);
  TString fo_name = TString::Format("../rootfile/data%04d_tw_doppler_%dns_jump_%dns.root", run, win, jump);
#else
  TString fi_name = TString::Format("../rootfile/data%04d_build_%dns_jump_%dns.root", run, win, jump);
  TString fo_name = TString::Format("../rootfile/data%04d_doppler_%dns_jump_%dns.root", run, win, jump);
#endif

  std::cout << "analysis " << fi_name << std::endl;

  doppler *dp = new doppler(fi_name.Data(), fo_name.Data(), CONFIG);

  dp->Process();
  dp->SaveFile();

  return 0;
}
