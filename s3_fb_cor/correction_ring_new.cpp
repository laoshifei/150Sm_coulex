#include <sys/stat.h>
#include <sys/types.h>

//
void correction_single(ofstream &fo, int run, int sector, int ring, double x0, double x1, double rob=0.6)
{
  gStyle->SetOptFit(1111);

  mkdir("fig", 0777);
  mkdir(TString::Format("./fig/%04d",run).Data(), 0777);

  TFile *fi = TFile::Open(TString::Format("../rootfile/si/data%04d_gr_fast.root", run).Data());
  if(!fi || fi->IsZombie()){
    cout << "can not open rootfile." << endl;
    return;
  }

  TCanvas *c2 = new TCanvas(TString::Format("c2_for_ring_sector%02d_ring%02d",sector,ring).Data(), "", 0, 0, 600, 400);
  TCanvas *c3 = new TCanvas(TString::Format("c3_for_ring_sector%02d_ring%02d",sector,ring).Data(), "", 0, 0, 600, 400);
  TGraph *gr = (TGraph*)gFile->Get(TString::Format("hh_ring%02d_sector%02d",ring,sector).Data());
  if(!gr){
    cout << "cannot find graph in rootfile." << endl;
    return;
  }

  c2->cd();
  gr->SetTitle(Form("for_ring_sector%02d_ring%02d",sector,ring));
  gr->GetXaxis()->SetTitle("ring energy [keV]");
  gr->GetYaxis()->SetTitle("sector energy [keV]");
  gr->Draw("ap");

  gr->Fit("pol1", TString::Format("ROB=%lf",rob).Data(), "same", x0, x1);
  TF1 *tf = (TF1*)gROOT->GetListOfFunctions()->FindObject("pol1");

  cout << sector << " " << ring << " " << tf->GetParameter(0) << " " << tf->GetParameter(1) << endl;

  // check
  c3->cd();
  double xx, yy;
  TGraph *gr_residual = new TGraph(gr->GetN());
  for(int i=0;i<gr->GetN();i++){
    xx = gr->GetPointX(i);
    yy = gr->GetPointY(i);
    xx = yy-tf->GetParameter(0)-tf->GetParameter(1)*xx;
    gr_residual->SetPoint(i, xx, yy);
  }

  gr_residual->GetXaxis()->SetRangeUser(-5000, 5000);
  gr_residual->SetTitle(Form("residual_between_sector%02d_ring%02d",sector,ring));
  gr_residual->GetXaxis()->SetTitle("residual error [keV]");
  gr_residual->GetYaxis()->SetTitle("sector energy [keV]");
  gr_residual->Draw("ap");

  c2->SaveAs(TString::Format("./fig/%04d/%s.png", run, c2->GetName()));
  c3->SaveAs(TString::Format("./fig/%04d/%s.png", run, c3->GetName()));

  fi->Close();

  fo << ring << "    " << tf->GetParameter(0) << "    " << tf->GetParameter(1) << endl;
}

//
void correction_ring_new(int run, double x0, double x1, double rob)
{
  gROOT->SetBatch(1);

  gSystem->Unlink(TString::Format("./par/correction_ring_%04d.txt",run).Data());
  ofstream fo(TString::Format("./par/correction_ring_%04d.txt",run).Data(), std::ios::app);
  fo << " ring  p0                p1" << endl;

  for(int i=1;i<=24;i++){
    correction_single(fo, run, 1, i, x0, x1, rob);
  }

  fo.close();
}
