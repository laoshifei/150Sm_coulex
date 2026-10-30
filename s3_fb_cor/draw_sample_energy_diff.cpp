//
void draw_sample_energy_diff(int run, int ring, int sector)
{
  TFile *fi = TFile::Open(Form("../rootfile/si/data%04d_build_s3_200ns.root",run));
  if(!fi || fi->IsZombie()) { Error("draw_sample_energy_diff", "Cannot open input"); return; }
  TTree *tr_event = (TTree*)fi->Get("tr_event");
  if(!tr_event) { Error("draw_sample_energy_diff", "Missing tr_event"); fi->Close(); return; }

  //
  TH1D *h_diff = new TH1D("h_diff", "", 5000, -25000, 25000);
  TString cut_str = Form("n_s3_sector==1 && n_s3_ring==1 && s3_ring_id==%d && s3_sector_id==%d", ring, sector);

  tr_event->Draw("s3_sector_energy-s3_ring_energy>>h_diff", cut_str.Data(), "goff", tr_event->GetEntries()/32);
  h_diff->SetTitle(Form("run %04d, ring %d / sector %d, first 1/32", run, ring, sector));
  h_diff->GetXaxis()->SetTitle("sector - ring energy [keV]");
  h_diff->GetYaxis()->SetTitle("counts");

  TCanvas *cc_diff = new TCanvas("cc_diff", "", 0, 0, 600, 400);
  cc_diff->cd();
  cc_diff->SetLogy();
  h_diff->Draw("");

  gSystem->mkdir("fig", true);
  cc_diff->SaveAs(Form("fig/run%04d_ring%02d_sector%02d_diff.png",run,ring,sector));
  fi->Close();
}
