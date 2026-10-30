//
void draw_sample_rs_th2(int run, int ring, int sector)
{
  //
  TFile *fi = TFile::Open(Form("../rootfile/si/data%04d_build_s3_200ns.root",run));
  if(!fi || fi->IsZombie()) { Error("draw_sample_rs_th2", "Cannot open old input"); return; }
  TTree *tr_event = (TTree*)fi->Get("tr_event");
  if(!tr_event) { Error("draw_sample_rs_th2", "Missing old tr_event"); fi->Close(); return; }

  TH2D *hh = new TH2D("hh", "", 1000, 0, 100000, 1000, 0, 100000);

  TString cut_str = Form("n_s3_ring==1 && n_s3_sector==1 && s3_ring_id==%d && s3_sector_id==%d", ring, sector);
  tr_event->Draw("s3_sector_energy:s3_ring_energy>>hh", cut_str.Data(), "goff", tr_event->GetEntries());
  hh->SetTitle("th2 between sector and ring");
  hh->GetXaxis()->SetTitle("ring energy [keV]");
  hh->GetYaxis()->SetTitle("sector energy [keV]");

  TCanvas *cc = new TCanvas("cc", "", 0, 0, 600, 400);
  cc->cd();
  cc->SetLogz();
  hh->Draw("colz");

  gSystem->mkdir("fig", true);
  cc->SaveAs(Form("fig/run%04d_ring%02d_sector%02d_previous.png",run,ring,sector));

  //
  TFile *fi_new = TFile::Open(Form("../rootfile/si/data%04d_build_s3_200ns_new.root",run));
  if(!fi_new || fi_new->IsZombie()) { Error("draw_sample_rs_th2", "Cannot open new input"); fi->Close(); return; }
  TTree *tr_event_new = (TTree*)fi_new->Get("tr_event");
  if(!tr_event_new) { Error("draw_sample_rs_th2", "Missing new tr_event"); fi->Close(); fi_new->Close(); return; }

  TH2D *hh_new = new TH2D("hh_new", "", 1000, 0, 100000, 1000, 0, 100000);

  tr_event_new->Draw("s3_sector_energy:s3_ring_energy>>hh_new", cut_str.Data(), "goff", tr_event_new->GetEntries());
  hh_new->SetTitle("th2 between sector and ring new");
  hh_new->GetXaxis()->SetTitle("ring energy [keV]");
  hh_new->GetYaxis()->SetTitle("sector energy [keV]");

  TCanvas *cc_new = new TCanvas("cc_new", "", 0, 0, 600, 400);
  cc_new->cd();
  cc_new->SetLogz();
  hh_new->Draw("colz");

  cc_new->SaveAs(Form("fig/run%04d_ring%02d_sector%02d_candidate.png",run,ring,sector));
  fi->Close(); fi_new->Close();
}
