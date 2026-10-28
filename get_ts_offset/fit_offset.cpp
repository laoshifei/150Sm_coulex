#include "offset_fit.h"
// root -l -b -q 'fit_offset.cpp("candidate.root",617,"candidate_fit","isolated",100,100.,"legacy")'
// scope="particle"：逐粒子通道对诊断，不自动解成或安装全局 offset。
int fit_offset(TString filename,int run,TString outdir="offset_fit_candidates",
               TString selection="isolated",int min_entries=100,double halfwidth_ns=100.,
               TString scope="legacy") {
  return sm150_offset::export_offsets(filename,run,outdir,true,selection.Data(),min_entries,halfwidth_ns,scope.Data());
}
