#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>

// Draw saved histograms only. No event processing or ROOT output mutation.
void draw_sample(const char* input,const char* output_directory) {
  namespace fs=std::filesystem;
  std::unique_ptr<TFile> file(TFile::Open(input,"READ"));
  if(!file||file->IsZombie())throw std::runtime_error("Cannot read analysis ROOT");
  const fs::path dir(output_directory);
  if(fs::exists(dir)||fs::is_symlink(dir))throw std::runtime_error("Choose a new figure directory");
  const char* names[]={"multiplicity_clean","face_de_energy","particle_dt","phi_after_time",
                      "cm_after_phi","energy_low_high_before","energy_low_high_after",
                      "energy_theta_before","energy_theta_after","energy_low_after","energy_high_after"};
  for(auto name:names)if(!dynamic_cast<TH1*>(file->Get(name)))throw std::runtime_error(std::string("Missing histogram: ")+name);
  fs::create_directories(dir);gROOT->SetBatch(true);gStyle->SetOptStat(0);
  for(auto name:names) {
    auto* h=dynamic_cast<TH1*>(file->Get(name));
    TCanvas canvas(("canvas_"+std::string(name)).c_str(),name,900,650);
    canvas.SetRightMargin(h->GetDimension()==2?.15:.05);
    if(h->GetMaximum()>0) {
      if(h->GetDimension()==2)canvas.SetLogz();
      else {canvas.SetLogy();h->SetMinimum(.5);}
    }
    h->Draw(h->GetDimension()==2?"colz":"hist");
    canvas.SaveAs((dir/(std::string(name)+".png")).c_str());
  }
}
