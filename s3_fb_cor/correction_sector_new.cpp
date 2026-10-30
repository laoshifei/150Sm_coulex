#include <sys/stat.h>
#include <sys/types.h>

//
void correction_single(ofstream &fo, int run, int sector, int ring, double x0, double x1, double rob=0.6)
{
  mkdir("fig", 0777);
  mkdir(TString::Format("./fig/%04d", run).Data(), 0777);

  //
  ifstream fi_ring_cor(TString::Format("./par/correction_ring_%04d.txt",run).Data());
  if(!fi_ring_cor){
    cout << "can not open ring correction file." << endl;
    return;
  }

  std::map<int, std::vector<double>> map_ring_cor;
  std::string line;
  std::getline(fi_ring_cor, line);

  int id;
  double p0, p1;
  while(1){
    fi_ring_cor >> id >> p0 >> p1;
    if(!fi_ring_cor.good()) break;

    map_ring_cor[id] = {p0, p1};
  }

  for(auto &it : map_ring_cor){
    cout << it.first << " " << it.second[0] << " " << it.second[1] << endl;
  }
  fi_ring_cor.close();

  //
  TFile *fi = TFile::Open(TString::Format("../rootfile/si/data%04d_gr_fast.root", run).Data());
  if(!fi || fi->IsZombie()){
    cout << "can not open rootfile." << endl;
    return;
  }

  TCanvas *c2 = new TCanvas(TString::Format("c2_for_sector_sector%02d_ring%02d",sector,ring).Data(), "", 0, 0, 600, 400);
  TCanvas *c3 = new TCanvas(TString::Format("c3_for_sector_sector%02d_ring%02d",sector,ring).Data(), "", 0, 0, 600, 400);
  
  TGraph *gr = (TGraph*)gFile->Get(TString::Format("hh_sector%02d_ring%02d",sector,ring).Data());
  if(!gr){
    cout << "cannot find graph in rootfile." << endl;
    return;
  }

  double xx, yy;
  TGraph *gr_new = new TGraph(gr->GetN());
  for(int i=0;i<gr->GetN();i++){
    xx = gr->GetPointX(i);
    yy = gr->GetPointY(i);
    yy = map_ring_cor[ring][0]+map_ring_cor[ring][1]*yy;
    gr_new->SetPoint(i, xx, yy);
  }

  c2->cd();
  gr_new->SetTitle(Form("for_sector_sector%02d_ring%02d",sector,ring));
  gr_new->GetXaxis()->SetTitle("sector energy [keV]");
  gr_new->GetYaxis()->SetTitle("corrected ring energy [keV]");
  gr_new->Draw("ap");

  gr_new->Fit("pol1", TString::Format("ROB=%lf",rob).Data(), "same", x0, x1);
  TF1 *tf = (TF1*)gROOT->GetListOfFunctions()->FindObject("pol1");

  cout << sector << " " << ring << " " << tf->GetParameter(0) << " " << tf->GetParameter(1) << endl;

  // check
  c3->cd();
  TGraph *gr_residual = new TGraph(gr_new->GetN());
  for(int i=0;i<gr_new->GetN();i++){
    xx = gr_new->GetPointX(i);
    yy = gr_new->GetPointY(i);
    xx = yy-tf->GetParameter(0)-tf->GetParameter(1)*xx;
    gr_residual->SetPoint(i, xx, yy);
  }

  gr_residual->GetXaxis()->SetRangeUser(-5000, 5000);
  gr_residual->SetTitle(Form("residual_between_sector%02d_ring%02d",sector,ring));
  gr_residual->GetXaxis()->SetTitle("residual error [keV]");
  gr_residual->GetYaxis()->SetTitle("corrected ring energy [keV]");
  gr_residual->Draw("ap");

  c2->SaveAs(TString::Format("./fig/%04d/%s.png", run, c2->GetName()));
  c3->SaveAs(TString::Format("./fig/%04d/%s.png", run, c3->GetName()));

  fi->Close();

  fo << sector << "    " << tf->GetParameter(0) << "    " << tf->GetParameter(1) << endl;
}

//
void correction_sector_new(int run, double x0, double x1, double rob)
{
  gROOT->SetBatch(1);

  gSystem->Unlink(TString::Format("./par/correction_sector_%04d.txt",run).Data());
  ofstream fo(TString::Format("./par/correction_sector_%04d.txt",run).Data(), std::ios::app);
  fo << " sector  p0                p1" << endl;

  for(int i=1;i<=32;i++){
    correction_single(fo, run, i, 11, x0, x1, rob);
  }

  fo.close();
}
