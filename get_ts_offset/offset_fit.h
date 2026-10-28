#ifndef SM150_OFFSET_FIT_H
#define SM150_OFFSET_FIT_H
// Share export and quality checks between the two ROOT macros; write candidates with explicit status, not active parameters.
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TObjString.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TString.h>
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
namespace sm150_offset {
using json=nlohmann::json;
namespace fs=std::filesystem;
inline json estimate(TH1D& source,bool fit,int min_entries,double halfwidth) {
  json r={{"status","insufficient_statistics"},{"offset_ns",nullptr},{"error_ns",nullptr},
    {"entries",source.GetEntries()},{"method",fit?"gaussian_core_plus_constant":"peak_bin_center"},
    {"bin_ns",source.GetBinWidth(1)},{"fit_status",nullptr},{"covariance_status",nullptr},
    {"root_chi2",nullptr},{"ndf",nullptr},{"sigma_ns",nullptr},{"fit_range_ns",nullptr}};
  if(source.GetEntries()<min_entries || source.Integral()<min_entries) return r;
  std::unique_ptr<TH1D> h(static_cast<TH1D*>(source.Clone((std::string(source.GetName())+"_fit").c_str())));
  h->SetDirectory(nullptr);
  const int n=h->GetNbinsX(),peak=h->GetMaximumBin();
  const double center=h->GetBinCenter(peak),height=h->GetBinContent(peak),bw=h->GetBinWidth(peak);
  if(peak<=2 || peak>=n-1) { r["status"]="peak_at_search_edge"; return r; }
  std::vector<double> background;
  for(int i=1;i<=n;++i) if(std::abs(h->GetBinCenter(i)-center)>halfwidth) background.push_back(h->GetBinContent(i));
  if(background.empty()) { r["status"]="no_background_region"; return r; }
  std::sort(background.begin(),background.end());
  const double bg=background[background.size()/2];
  if(height-bg<5*std::sqrt(std::max(1.,bg))) { r["status"]="no_clear_peak"; return r; }
  if(!fit) {
    r["status"]="peak_bin_estimate"; r["offset_ns"]=center;
    r["bin_quantization_ns"]=bw/std::sqrt(12.);
    r["error_note"]="bin quantization is not a statistical peak-center uncertainty";
    return r;
  }
  double lo=std::max(h->GetXaxis()->GetXmin(),center-halfwidth);
  double hi=std::min(h->GetXaxis()->GetXmax(),center+halfwidth);
  if(hi-lo<12*bw) { r["status"]="insufficient_fit_bins"; return r; }
  TF1 function((std::string(source.GetName())+"_core").c_str(),"gaus(0)+pol0(3)",lo,hi);
  double sigma=std::max(2*bw,halfwidth/6);
  function.SetParameters(height-bg,center,sigma,bg);
  function.SetParLimits(0,0.,std::max(10.,height*20));
  function.SetParLimits(1,lo+2*bw,hi-2*bw);
  function.SetParLimits(2,bw/2,halfwidth/2);
  function.SetParLimits(3,0.,std::max(10.,height*2));
  // Poisson likelihood on original nonnegative counts, integrated bin model.
  TFitResultPtr result=h->Fit(&function,"QRS0LNI");
  const int status=int(result),cov=result.Get()?result->CovMatrixStatus():-1;
  const double mu=function.GetParameter(1),err=function.GetParError(1),s=function.GetParameter(2);
  r["fit_status"]=status;r["covariance_status"]=cov;r["fit_range_ns"]={lo,hi};
  r["root_chi2"]=function.GetChisquare();r["ndf"]=function.GetNDF();r["sigma_ns"]=s;
  r["fit_parameters"]={function.GetParameter(0),mu,s,function.GetParameter(3)};
  r["fitted_offset_ns"]=std::isfinite(mu)?json(mu):json(nullptr);
  r["fitted_error_ns"]=std::isfinite(err)?json(err):json(nullptr);
  bool valid=status==0 && cov==3 && function.GetNDF()>0 && std::isfinite(mu) &&
             std::isfinite(err) && err>0 && s>bw*.501 && s<halfwidth*.499 &&
             mu>lo+2.01*bw && mu<hi-2.01*bw;
  if(!valid) { r["status"]="fit_failed_or_boundary"; return r; }
  if(function.GetChisquare()/function.GetNDF()>5.) { r["status"]="core_shape_mismatch"; return r; }
  for(int i=2;i<n;++i)
    if(std::abs(h->GetBinCenter(i)-mu)>halfwidth && h->GetBinContent(i)-bg>.8*(height-bg) &&
       h->GetBinContent(i)>=h->GetBinContent(i-1) && h->GetBinContent(i)>=h->GetBinContent(i+1)) {
      r["status"]="ambiguous_multiple_peaks"; return r;
    }
  r["status"]="valid_candidate";r["offset_ns"]=mu;r["error_ns"]=err;
  return r;
}
inline int export_offsets(const TString& filename,int run,const TString& destination,bool fit,
                          const std::string& selection="isolated",int min_entries=100,
                          double halfwidth=100.,const std::string& scope="legacy") {
  if(selection!="all" && selection!="isolated") throw std::runtime_error("Selection must be all or isolated");
  if(scope!="legacy" && scope!="particle") throw std::runtime_error("Scope must be legacy or particle");
  if(min_entries<1 || !std::isfinite(halfwidth) || halfwidth<=0) throw std::runtime_error("Invalid fit settings");
  std::unique_ptr<TFile> f(TFile::Open(filename.Data(),"READ"));
  if(!f || f->IsZombie()) throw std::runtime_error("Cannot open timing ROOT");
  auto* info=dynamic_cast<TObjString*>(f->Get("metadata_json"));
  if(!info) throw std::runtime_error("Missing v2 metadata; regenerate histograms with repaired analysis");
  json meta=json::parse(info->GetString().Data());
  if(meta.at("schema")!="sm150_offset_v2" || meta.at("status")!="complete" ||
     std::stoi(meta.at("run").get<std::string>())!=run) throw std::runtime_error("Incomplete/wrong-run timing ROOT");
  if(scope=="legacy" && !meta.at("reference_present").get<bool>())
    throw std::runtime_error("Reference channel has no selected normal-CFD hits; refusing zero-offset table");
  fs::path out(destination.Data());
  if(!out.parent_path().empty()) fs::create_directories(out.parent_path());
  if(!fs::create_directory(out)) throw std::runtime_error("Export directory exists; use a fresh candidate directory");
  gROOT->SetBatch(kTRUE);
  json rows=json::array(),counts=json::object();
  std::ofstream csv(out/"offset_candidates.csv");
  if(!csv) throw std::runtime_error("Cannot create candidate table");
  csv<<"histogram,reference_key,target_key,category,mode,selection,status,offset_ns,error_ns,entries,bin_ns,fit_status,covariance_status,root_chi2,ndf,sigma_ns\n";
  for(auto it=meta.at("histograms").begin();it!=meta.at("histograms").end();++it) {
    const std::string category=it.value().at("category");
    const bool legacy=category=="ge_ge" || category=="particle_ge";
    if((scope=="legacy")!=legacy || category=="ge_particle") continue;
    std::string id=it.key(),name=id+(selection=="isolated"?"__isolated":"");
    auto* h=dynamic_cast<TH1D*>(f->Get(name.c_str()));
    if(!h) throw std::runtime_error("Missing histogram "+name);
    const int ref=it.value().at("reference_key"),key=it.value().at("target_key");
    json r=estimate(*h,fit,min_entries,halfwidth);
    if(legacy && key==ref) { r["status"]="reference";r["offset_ns"]=0.;r["error_ns"]=nullptr; }
    r["histogram"]=id;r["reference_key"]=ref;r["target_key"]=key;r["category"]=category;
    r["mode"]=meta.at("mode");r["selection"]=selection;
    std::string status=r.at("status"); counts[status]=counts.value(status,0)+1;
    const std::vector<std::string> columns={"histogram","reference_key","target_key","category","mode","selection",
      "status","offset_ns","error_ns","entries","bin_ns","fit_status","covariance_status","root_chi2","ndf","sigma_ns"};
    for(size_t i=0;i<columns.size();++i) {
      if(i) csv<<',';
      const auto& v=r.at(columns[i]); if(v.is_string()) csv<<v.get<std::string>(); else if(!v.is_null()) csv<<v.dump();
    }
    csv<<'\n'; rows.push_back(r);
    if(h->GetEntries()>=min_entries) {
      TCanvas canvas(("canvas_"+id).c_str(),id.c_str(),1000,500); h->Draw("hist");
      // Retain the full time range; the red curve describes only the fitted core and does not imply that the broad tail is corrected.
      std::unique_ptr<TF1> curve;
      if(r.contains("fit_parameters")) {
        auto v=r["fit_parameters"].get<std::vector<double>>();
        curve=std::make_unique<TF1>((id+"_draw").c_str(),"gaus(0)+pol0(3)",
                                  r["fit_range_ns"][0].get<double>(),r["fit_range_ns"][1].get<double>());
        curve->SetParameters(v.data());curve->SetLineColor(kRed);curve->Draw("same");
      }
      canvas.SaveAs((out/(id+".png")).c_str());
    }
  }
  csv.close(); if(!csv) throw std::runtime_error("Candidate CSV write failed");
  json report={{"status","complete"},{"source",fs::absolute(filename.Data()).string()},
    {"run",run},{"mode",meta.at("mode")},{"selection",selection},{"scope",scope},{"min_entries",min_entries},
    {"fit_halfwidth_ns",halfwidth},{"status_counts",counts},{"rows",rows},
    {"offset_interpretation","dt=reference-target: add fitted center to target relative to fixed reference; residual is a remaining term, not a new total"},
    {"adopted",false},{"conditions","Core Gaussian+constant only; thresholds are screening choices. Combined CFD spectra can remain non-Gaussian."}};
  std::ofstream js(out/"offset_candidates.json");js<<report.dump(2)<<'\n';js.close();
  if(!js) throw std::runtime_error("Candidate JSON write failed");
  return 0;
}
} // namespace
#endif
