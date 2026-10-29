#include "sort.h"

#include <iostream>
#include <fstream>
#include <assert.h>
#include <stdlib.h>

#include "TString.h"
#include "TStopwatch.h"

//
int main(int argc, char const *argv[])
{
  if(argc != 2 && argc != 4){
    std::cout << "./sort run [input.root output.root]" << std::endl;
    return -1;
  }

  int run = atoi(argv[1]);
  TString file_in = TString::Format("../rootfile/data_C2_%04d_%04d.root", run, run);

#ifdef TW
  TString file_out = TString::Format("../rootfile/data%04d_tw.root", run);
#else
  TString file_out = TString::Format("../rootfile/data%04d.root", run);
#endif
  if(argc == 4){
    file_in = argv[2];
    file_out = argv[3];
  }
  if(file_in == file_out){
    std::cerr << "Input and output ROOT files must be different." << std::endl;
    return -1;
  }
  std::cout << "input " << file_in << " -> output " << file_out << std::endl;

  TStopwatch timer;
  timer.Start();

  sort *so = new sort(file_in.Data(), file_out.Data(), run);
  so->Process();

  timer.Stop();
  Double_t real_time = timer.RealTime();
  Double_t cpu_time = timer.CpuTime();

  printf("time : RT=%7.3f s, Cpu=%7.3f s\n", real_time, cpu_time);

  delete so;

  return 0;
}
