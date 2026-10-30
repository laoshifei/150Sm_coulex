// Redraw small histogram outputs; no decoded event input or recalibration.
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TKey.h>
#include <TObjString.h>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>

void draw_sample(const char* input,const char* new_output_dir) {
  namespace fs=std::filesystem;
  const fs::path output=fs::absolute(new_output_dir);
  if(fs::exists(output)||fs::is_symlink(output))throw std::runtime_error("Figure directory already exists");
  std::unique_ptr<TFile> f(TFile::Open(input,"READ"));
  if(!f||f->IsZombie())throw std::runtime_error("Cannot open histogram file");
  auto* meta=dynamic_cast<TObjString*>(f->Get("metadata"));
  if(!meta||!meta->GetString().Contains("sm150_raw2hist_v1"))throw std::runtime_error("Not a current raw2hist output");
  fs::create_directories(output);TCanvas canvas("raw2hist_channel","Channel singles",1000,750);
  TIter next(f->GetListOfKeys());
  while(auto* key=dynamic_cast<TKey*>(next())) {
    auto* dir=f->GetDirectory(key->GetName());if(!dir)continue;
    canvas.Clear();canvas.Divide(1,2);
    auto* adc=dynamic_cast<TH1D*>(dir->Get("adc_clean"));
    auto* energy=dynamic_cast<TH1D*>(dir->Get("energy_final"));
    canvas.cd(1);gPad->SetLogy();if(adc){adc->SetTitle((std::string(key->GetName())+";ADC channel;Counts / 2 ADC channels").c_str());adc->SetMinimum(.7);adc->Draw("hist");}
    canvas.cd(2);gPad->SetLogy();if(energy){energy->SetTitle((std::string(key->GetName())+";Energy (keV);Counts / bin").c_str());energy->SetMinimum(.7);energy->Draw("hist");}
    canvas.SaveAs((output/(std::string(key->GetName())+".png")).c_str());
  }
}
