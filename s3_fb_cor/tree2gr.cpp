//
void tree2gr_sector(int run)
{
  TFile *fi = TFile::Open(TString::Format("../rootfile/si/data%04d_build_s3_200ns_no_s3cor.root", run).Data());
  if(!fi || fi->IsZombie()){
    cout << "can not open rootfile." << endl;
    return;
  }

  TTree *tr = (TTree*)fi->Get("tr_event");
  if(!tr) { Error("legacy_s3", "Missing tr_event"); fi->Close(); return; }
  TH2D *hh_for_sector[32];
  TGraph *gr[32];

  int n = 0;
  TString draw_str, cut_str;
  for(int i=1;i<=32;i++){
    hh_for_sector[i-1] = new TH2D(TString::Format("hh_sector%02d_ring%02d",i,11).Data(), "", 8192, 0, 65536, 8192, 0, 65536);
    TString draw_str = TString::Format("s3_sector_energy:s3_ring_energy>>%s", hh_for_sector[i-1]->GetName());
    TString cut_str = TString::Format("n_s3_sector==1 && n_s3_ring==1 && s3_sector_id==%d && s3_ring_id==11",i);
    cout << draw_str << endl;
    cout << cut_str << endl;
    cout << endl;

    n = tr->Draw(draw_str.Data(), cut_str.Data(), "goff");
    cout << "n " << n << endl;
    gr[i-1] = new TGraph(n, tr->GetV1(), tr->GetV2());
    gr[i-1]->SetName(hh_for_sector[i-1]->GetName());

    if(hh_for_sector[i-1]) delete hh_for_sector[i-1];
  }

  TFile *fo = new TFile(TString::Format("../rootfile/si/data%04d_gr2.root", run).Data(), "recreate");

  fo->cd();
  for(int i=0;i<32;i++){
    gr[i]->Write();
  }

  fo->Close();

  fi->Close();
}

//
void tree2gr_ring(int run)
{
  TFile *fi = TFile::Open(TString::Format("../rootfile/si/data%04d_build_s3_200ns_no_s3cor.root", run).Data());
  if(!fi || fi->IsZombie()){
    cout << "can not open rootfile." << endl;
    return;
  }

  TTree *tr = (TTree*)fi->Get("tr_event");
  if(!tr) { Error("legacy_s3", "Missing tr_event"); fi->Close(); return; }
  TH2D *hh_for_ring[24];
  TGraph *gr[24];

  int n = 0;
  TString draw_str, cut_str;
  for(int i=1;i<=24;i++){
    hh_for_ring[i-1] = new TH2D(TString::Format("hh_ring%02d_sector01",i).Data(), "", 8192, 0, 65536, 8192, 0, 65536);
    TString draw_str = TString::Format("s3_sector_energy:s3_ring_energy>>%s", hh_for_ring[i-1]->GetName());
    TString cut_str = TString::Format("n_s3_sector==1 && n_s3_ring==1 && s3_sector_id==1 && s3_ring_id==%d",i);
    cout << draw_str << endl;
    cout << cut_str << endl;
    cout << endl;

    n = tr->Draw(draw_str.Data(), cut_str.Data(), "goff");
    cout << "n " << n << endl;
    gr[i-1] = new TGraph(n, tr->GetV2(), tr->GetV1());
    gr[i-1]->SetName(hh_for_ring[i-1]->GetName());

    if(hh_for_ring[i-1]) delete hh_for_ring[i-1];
  }

  TFile *fo = new TFile(TString::Format("../rootfile/si/data%04d_gr2.root", run).Data(), "update");

  fo->cd();
  for(int i=0;i<24;i++){
    gr[i]->Write();
  }

  fo->Close();

  fi->Close();
}

//
void tree2gr(int run)
{
  tree2gr_sector(run);
  tree2gr_ring(run);
}
