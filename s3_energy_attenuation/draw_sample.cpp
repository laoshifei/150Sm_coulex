// View native extraction output; no fitting or event processing.
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TROOT.h>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <memory>
#include <stdexcept>
void draw_sample(const char* input,const char* output) {
  namespace fs=std::filesystem;
  gROOT->SetBatch(true);
  if(fs::exists(output))throw std::runtime_error("Figure directory already exists");
  std::unique_ptr<TFile> source(TFile::Open(input,"READ"));
  if(!source || source->IsZombie())throw std::runtime_error("Cannot read spectra");
  auto* meta=dynamic_cast<TNamed*>(source->Get("metadata"));
  if(!meta)throw std::runtime_error("Missing source metadata");
  auto info=nlohmann::json::parse(meta->GetTitle());
  if(info.at("schema")!="sm150_s3_attenuation_spectra_v1")throw std::runtime_error("Wrong input schema");
  fs::create_directories(output);
  TCanvas canvas("s3_energy","S3 singles before attenuation",900,600);canvas.SetLogy();
  for(const auto& channel:info.at("channels")) {
    const std::string name=channel.at("name"),run=info.at("run");
    auto* h=dynamic_cast<TH2D*>(source->Get(("clean_"+name).c_str()));
    if(!h)throw std::runtime_error("Missing channel");
    std::unique_ptr<TH1D> projection(h->ProjectionX(("view_"+name).c_str(),1,h->GetNbinsY()));
    projection->SetDirectory(nullptr);
    projection->SetTitle(("Run "+run+" / "+name+";Energy before attenuation [keV];Clean singles / 50 keV").c_str());
    projection->SetMinimum(.5);projection->Draw("HIST");
    canvas.SaveAs((fs::path(output)/(name+".png")).c_str());
  }
}
