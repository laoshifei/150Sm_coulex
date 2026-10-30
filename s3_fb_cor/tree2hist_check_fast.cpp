//
void tree2hist_check_fast(int run)
{
  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);

  TFile *fi = TFile::Open(TString::Format("../rootfile/si/data%04d_build_s3_200ns.root", run).Data());
  if(!fi || fi->IsZombie()){
    cout << "can not open rootfile." << endl;
    return;
  }

  TTree *tr = (TTree*)fi->Get("tr_event");
  if(!tr) { Error("legacy_s3", "Missing tr_event"); fi->Close(); return; }

  //
  TH2D *hh_sector = new TH2D("hh_sector", "sector vs. id TH2", 32, 1, 33, 1000, 0, 100000);
  tr->Draw("s3_sector_energy:s3_sector_id>>hh_sector", "n_s3_sector==1 && n_s3_ring==1", "goff", tr->GetEntries()/32);
  TCanvas *cc_s = new TCanvas("cc_s", "", 0, 0, 600, 400);
  cc_s->cd();
  cc_s->SetLogz();
  hh_sector->GetXaxis()->SetTitle("id");
  hh_sector->GetYaxis()->SetTitle("energy [keV]");
  hh_sector->Draw("colz");
  cc_s->SaveAs(Form("fig/%04d/hh_sector.png", run));

  //
  TH2D *hh_ring = new TH2D("hh_ring", "ring vs. id TH2", 24, 1, 25, 1000, 0, 100000);
  tr->Draw("s3_ring_energy:s3_ring_id>>hh_ring", "n_s3_ring==1 && n_s3_sector==1", "goff", tr->GetEntries()/32);
  TCanvas *cc_r = new TCanvas("cc_r", "", 0, 0, 600, 400);
  cc_r->cd();
  cc_r->SetLogz();
  hh_ring->GetXaxis()->SetTitle("id");
  hh_ring->GetYaxis()->SetTitle("energy [keV]");
  hh_ring->Draw("colz");
  cc_r->SaveAs(Form("fig/%04d/hh_ring.png", run));

  //
  TH1D *h_sector[32], *h_sector2[32], *h_ring[24], *h_ring2[24];

  int n = 0;
  TCanvas *cc_sector[32];
  TLegend *leg_sector[32];
  for(int i=0;i<32;i++){
    cc_sector[i] = new TCanvas(Form("cc_sector%02d",i+1), "", 0, 0, 600, 400);

    //
    h_sector[i] = new TH1D(Form("h_sector%02d",i+1), "", 1000, 0, 100000);

    TString draw_str_sector = TString::Format("s3_sector_energy>>%s", h_sector[i]->GetName());
    TString cut_str_sector = TString::Format("n_s3_sector==1 && n_s3_ring==1 && s3_sector_id==%d",i+1);
    cout << draw_str_sector << endl;
    cout << cut_str_sector << endl;
    cout << endl;

    n = tr->Draw(draw_str_sector.Data(), cut_str_sector.Data(), "goff", tr->GetEntries()/32);
    cout << "n " << n << endl;

    //
    h_sector2[i] = new TH1D(Form("h_sector%02d_coin_ring",i+1), "", 1000, 0, 100000);

    TString draw_str_sector2 = TString::Format("s3_ring_energy>>%s", h_sector2[i]->GetName());
    TString cut_str_sector2 = TString::Format("n_s3_sector==1 && n_s3_ring==1 && s3_sector_id==%d",i+1);
    cout << draw_str_sector2 << endl;
    cout << cut_str_sector2 << endl;
    cout << endl;

    n = tr->Draw(draw_str_sector2.Data(), cut_str_sector2.Data(), "goff", tr->GetEntries()/32);
    cout << "n " << n << endl;

    //
    cc_sector[i]->cd();
    cc_sector[i]->SetLogy();
    h_sector[i]->SetLineColor(2);
    h_sector2[i]->SetLineColor(4);

    h_sector[i]->SetTitle(Form("sector%02d coincidence", i+1));
    h_sector[i]->GetXaxis()->SetTitle("energy [keV]");
    h_sector[i]->GetYaxis()->SetTitle("counts");
    h_sector[i]->Draw();
    h_sector2[i]->Draw("same");

    //
    leg_sector[i] = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg_sector[i]->AddEntry(h_sector[i], "sector", "l");
    leg_sector[i]->AddEntry(h_sector2[i], "ring", "l");
    leg_sector[i]->Draw("same");

    cc_sector[i]->SaveAs(Form("fig/%04d/%s.png", run, cc_sector[i]->GetName()));
  }

  TCanvas *cc_ring[24];
  TLegend *leg_ring[24];
  for(int i=0;i<24;i++){
    cc_ring[i] = new TCanvas(Form("cc_ring%02d",i+1), "", 0, 0, 600, 400);

    //
    h_ring[i] = new TH1D(Form("h_ring%02d",i+1), "", 1000, 0, 100000);

    TString draw_str_ring = TString::Format("s3_ring_energy>>%s", h_ring[i]->GetName());
    TString cut_str_ring = TString::Format("n_s3_ring==1 && n_s3_sector==1 && s3_ring_id==%d",i+1);
    cout << draw_str_ring << endl;
    cout << cut_str_ring << endl;
    cout << endl;

    n = tr->Draw(draw_str_ring.Data(), cut_str_ring.Data(), "goff", tr->GetEntries()/32);
    cout << "n " << n << endl;

    //
    h_ring2[i] = new TH1D(Form("h_ring%02d_coin_ring",i+1), "", 1000, 0, 100000);

    TString draw_str_ring2 = TString::Format("s3_sector_energy>>%s", h_ring2[i]->GetName());
    TString cut_str_ring2 = TString::Format("n_s3_ring==1 && n_s3_sector==1 && s3_ring_id==%d",i+1);
    cout << draw_str_ring2 << endl;
    cout << cut_str_ring2 << endl;
    cout << endl;

    n = tr->Draw(draw_str_ring2.Data(), cut_str_ring2.Data(), "goff", tr->GetEntries()/32);
    cout << "n " << n << endl;

    //
    cc_ring[i]->cd();
    cc_ring[i]->SetLogy();
    h_ring[i]->SetLineColor(2);
    h_ring2[i]->SetLineColor(4);

    h_ring[i]->SetTitle(Form("ring%02d coincidence", i+1));
    h_ring[i]->GetXaxis()->SetTitle("energy [keV]");
    h_ring[i]->GetYaxis()->SetTitle("counts");
    h_ring[i]->Draw();
    h_ring2[i]->Draw("same");

    //
    leg_ring[i] = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg_ring[i]->AddEntry(h_ring[i], "ring", "l");
    leg_ring[i]->AddEntry(h_ring2[i], "sector", "l");
    leg_ring[i]->Draw("same");

    cc_ring[i]->SaveAs(Form("fig/%04d/%s.png", run, cc_ring[i]->GetName()));
  }

  fi->Close();
}


