#include <sys/stat.h>
#include <sys/types.h>

//
void draw_sample_single()
{
  TFile *fi = TFile::Open("../rootfile/si/data0616_build_s3_200ns_no_s3cor.root");
  if(!fi || fi->IsZombie()) { Error("draw_sample_single", "Cannot open input"); return; }
  TTree *tr_event = (TTree*)fi->Get("tr_event");
  if(!tr_event) { Error("draw_sample_single", "Missing tr_event"); fi->Close(); return; }

  TH2D *hh = new TH2D("hh", "", 500,0,100000,500,0,100000);
  hh->GetXaxis()->SetTitle("ring energy [keV]");
  hh->GetYaxis()->SetTitle("sector energy [keV]");

  Long64_t n = tr_event->Draw("s3_sector_energy:s3_ring_energy>>hh", "", "goff", tr_event->GetEntries()/32);
  cout << "n " << n << endl;

  TCanvas *cc = new TCanvas("cc", "", 0, 0, 600, 400);
  cc->cd();
  cc->SetLogz();
  hh->Draw("colz");
  gSystem->mkdir("fig", true);
  cc->SaveAs("fig/cc.png");
  fi->Close();
}

//
void draw_sector(int run)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //
  gSystem->mkdir(Form("./fig/%04d",run), true);

  //
  TFile *fi;
  TGraph *hh[32];
  fi = TFile::Open(TString::Format("../rootfile/si/data%04d_gr_fast.root", run).Data());
  if(!fi || fi->IsZombie()){
    cout << "wrong rootfile." << endl;
    return;
  }

  for(int i=0;i<32;i++){
    hh[i] = (TGraph*)gFile->Get(Form("hh_sector%02d_ring11", i+1));
  }

  TCanvas *cc[32];
  for(int i=0;i<32;i++){
    cc[i] = new TCanvas(Form("cc_sector%02d_vs_ring11",i+1), Form("cc_ring%02d_vs_sector01",i+1), 0, 0, 600, 400);
    cc[i]->cd();
    cc[i]->SetLeftMargin(0.15);
    
    if(!hh[i]) { Error("draw_sample", "Missing channel graph"); continue; }
    hh[i]->GetXaxis()->SetLimits(0, 100000);
    hh[i]->SetMinimum(0); hh[i]->SetMaximum(100000);
    hh[i]->GetXaxis()->SetTitle(Form("sector%02d [keV]",i+1));
    hh[i]->GetYaxis()->SetTitle("ring11 [keV]");
    hh[i]->Draw("ap");  

    cc[i]->SaveAs(Form("fig/%04d/%s.png",run,cc[i]->GetName()));
  }

  fi->Close();
}

//
void draw_ring(int run)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //
  gSystem->mkdir(Form("./fig/%04d",run), true);

  //
  TFile *fi;
  TGraph *hh[24];
  fi = TFile::Open(TString::Format("../rootfile/si/data%04d_gr_fast.root", run).Data());
  if(!fi || fi->IsZombie()){
    cout << "wrong rootfile." << endl;
    return;
  }

  for(int i=0;i<24;i++){
    hh[i] = (TGraph*)gFile->Get(Form("hh_ring%02d_sector01", i+1));
  }

  TCanvas *cc[24];
  for(int i=0;i<24;i++){
    cc[i] = new TCanvas(Form("cc_ring%02d_vs_sector01",i+1), Form("cc_ring%02d_vs_sector01",i+1), 0, 0, 600, 400);
    cc[i]->cd();
    cc[i]->SetLeftMargin(0.15);
    
    if(!hh[i]) { Error("draw_sample", "Missing channel graph"); continue; }
    hh[i]->GetXaxis()->SetLimits(0, 100000);
    hh[i]->SetMinimum(0); hh[i]->SetMaximum(100000);
    hh[i]->GetXaxis()->SetTitle(Form("ring%02d [keV]",i+1));
    hh[i]->GetYaxis()->SetTitle("sector1 [keV]");
    hh[i]->Draw("ap");  

    cc[i]->SaveAs(Form("fig/%04d/%s.png",run,cc[i]->GetName()));
  }

  fi->Close();
}

//
void draw_sample(int run)
{
  draw_sector(run);
  draw_ring(run);
}
