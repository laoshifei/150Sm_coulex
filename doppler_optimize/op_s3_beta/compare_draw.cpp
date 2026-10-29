void compare_draw()
{
  gStyle->SetOptStat(0);

  int colors[12] = {
    TColor::GetColor("#E41A1C"), TColor::GetColor("#377EB8"),
    TColor::GetColor("#4DAF4A"), TColor::GetColor("#984EA3"),
    TColor::GetColor("#FF7F00"), TColor::GetColor("#A65628"),
    TColor::GetColor("#F781BF"), TColor::GetColor("#999999"),
    TColor::GetColor("#66C2A5"), TColor::GetColor("#FC8D62"),
    TColor::GetColor("#8DA0CB"), TColor::GetColor("#FFD92F")
  };

  vector<string> v_filename = {
    "try/014/rootfile/data0616_0674_doppler_200ns_jump_300ns_hist.root",
    "try/012/rootfile/data0616_0674_doppler_200ns_jump_300ns_hist.root",
    "try/010/rootfile/data0616_0674_doppler_200ns_jump_300ns_hist.root",
    "try/008/rootfile/data0616_0674_doppler_200ns_jump_300ns_hist.root",
    "try/006/rootfile/data0616_0674_doppler_200ns_jump_300ns_hist.root",
    "try/004/rootfile/data0616_0674_doppler_200ns_jump_300ns_hist.root",
    "try/002/rootfile/data0616_0674_doppler_200ns_jump_300ns_hist.root",
    "try/001/rootfile/data0616_0674_doppler_200ns_jump_300ns_hist.root"
  };

  vector<TFile*> v_fi;
  vector<TH1D*> v_h_ga_s3a, v_hh_ga_s3a, v_h_grs32_s3a, v_hh_grs32_s3a, v_h_grs41_s3a, v_hh_grs41_s3a, v_h_grs42_s3a, v_hh_grs42_s3a, v_h_grs43_s3a, v_hh_grs43_s3a, v_h_grs52_s3a, v_hh_grs52_s3a, v_h_grs53_s3a, v_hh_grs53_s3a;

  for(int i=0;i<v_filename.size();++i){
    v_fi.push_back(TFile::Open(Form("%s", v_filename[i].c_str())));
    if(v_fi[i]->IsZombie()){
      cout << "wrong " << v_filename[i] << endl;
      return;
    }

    // s3 all
    v_h_ga_s3a.push_back((TH1D*)v_fi[i]->Get("event_e_dc_r_ge_all_s3_all"));
    v_hh_ga_s3a.push_back((TH1D*)v_fi[i]->Get("bg_e_dc_r_ge_all_s3_all"));
    v_h_ga_s3a[i]->Add(v_hh_ga_s3a[i], -1);

    v_h_grs32_s3a.push_back((TH1D*)v_fi[i]->Get("event_e_dc_r_ge_ring3_sector2_s3_all"));
    v_hh_grs32_s3a.push_back((TH1D*)v_fi[i]->Get("bg_e_dc_r_ge_ring3_sector2_s3_all"));
    v_h_grs32_s3a[i]->Add(v_hh_grs32_s3a[i], -1);

    v_h_grs41_s3a.push_back((TH1D*)v_fi[i]->Get("event_e_dc_r_ge_ring4_sector1_s3_all"));
    v_hh_grs41_s3a.push_back((TH1D*)v_fi[i]->Get("bg_e_dc_r_ge_ring4_sector1_s3_all"));
    v_h_grs41_s3a[i]->Add(v_hh_grs41_s3a[i], -1);

    v_h_grs42_s3a.push_back((TH1D*)v_fi[i]->Get("event_e_dc_r_ge_ring4_sector2_s3_all"));
    v_hh_grs42_s3a.push_back((TH1D*)v_fi[i]->Get("bg_e_dc_r_ge_ring4_sector2_s3_all"));
    v_h_grs42_s3a[i]->Add(v_hh_grs42_s3a[i], -1);

    v_h_grs43_s3a.push_back((TH1D*)v_fi[i]->Get("event_e_dc_r_ge_ring4_sector3_s3_all"));
    v_hh_grs43_s3a.push_back((TH1D*)v_fi[i]->Get("bg_e_dc_r_ge_ring4_sector3_s3_all"));
    v_h_grs43_s3a[i]->Add(v_hh_grs43_s3a[i], -1);

    v_h_grs52_s3a.push_back((TH1D*)v_fi[i]->Get("event_e_dc_r_ge_ring5_sector2_s3_all"));
    v_hh_grs52_s3a.push_back((TH1D*)v_fi[i]->Get("bg_e_dc_r_ge_ring5_sector2_s3_all"));
    v_h_grs52_s3a[i]->Add(v_hh_grs52_s3a[i], -1);

    v_h_grs53_s3a.push_back((TH1D*)v_fi[i]->Get("event_e_dc_r_ge_ring5_sector3_s3_all"));
    v_hh_grs53_s3a.push_back((TH1D*)v_fi[i]->Get("bg_e_dc_r_ge_ring5_sector3_s3_all"));
    v_h_grs53_s3a[i]->Add(v_hh_grs53_s3a[i], -1);
  }

  //
  int n_files = v_filename.size();

  TLine *l_334 = new TLine(334, 1., 334, v_h_ga_s3a[0]->GetMaximum()*10);
  l_334->SetLineColor(kRed);
  l_334->SetLineWidth(2);
  l_334->SetLineStyle(2);

  TLine *l_406 = new TLine(406.5, 1., 406.5, v_h_ga_s3a[0]->GetMaximum()*10);
  l_406->SetLineColor(kRed);
  l_406->SetLineWidth(2);
  l_406->SetLineStyle(2);

  TLine *l_439 = new TLine(439.4, 1., 439.4, v_h_ga_s3a[0]->GetMaximum()*10);
  l_439->SetLineColor(kRed);
  l_439->SetLineWidth(2);
  l_439->SetLineStyle(2);

  // ge all, s3 all
  TCanvas *cc_ga_s3a = new TCanvas("cc_ga_s3a","",900,500);
  cc_ga_s3a->cd();
  cc_ga_s3a->SetLogy();

  TLegend *leg_ga_s3a = new TLegend(0.75, 0.38, 0.88, 0.92);
  for(int i=0;i<n_files;++i){
    leg_ga_s3a->AddEntry(v_h_ga_s3a[i], v_filename[i].substr(4,3).c_str());
    v_h_ga_s3a[i]->SetLineColor(colors[i]);
    v_h_ga_s3a[i]->SetLineWidth(2);
    v_h_ga_s3a[i]->GetXaxis()->SetRangeUser(300,500);
    if(i==0){
      v_h_ga_s3a[i]->SetTitle("ge all, different beta pars");
      v_h_ga_s3a[i]->GetXaxis()->SetTitle("energy [keV]");
      v_h_ga_s3a[i]->GetYaxis()->SetTitle("counts / keV");
      v_h_ga_s3a[i]->GetYaxis()->SetRangeUser(1., v_h_ga_s3a[i]->GetMaximum()*10);
      v_h_ga_s3a[i]->Draw("hist");
    }else v_h_ga_s3a[i]->Draw("same hist");
  }
  leg_ga_s3a->Draw("same");

  l_334->Draw("same");
  l_406->Draw("same");
  l_439->Draw("same");

  cc_ga_s3a->SaveAs("fig/cc_ga_s3a.png");

  //
  // ge ring3 sector2, s3 all
  TCanvas *cc_grs32_s3a = new TCanvas("cc_grs32_s3a","",900,500);
  cc_grs32_s3a->cd();
  cc_grs32_s3a->SetLogy();

  TLegend *leg_grs32_s3a = new TLegend(0.75, 0.38, 0.88, 0.92);
  for(int i=0;i<n_files;++i){
    leg_grs32_s3a->AddEntry(v_h_grs32_s3a[i], v_filename[i].substr(4,3).c_str());
    v_h_grs32_s3a[i]->SetTitle("ge ring3 sector2, different beta pars");
    v_h_grs32_s3a[i]->SetLineColor(colors[i]);
    v_h_grs32_s3a[i]->SetLineWidth(2);
    v_h_grs32_s3a[i]->GetXaxis()->SetRangeUser(300,500);
    if(i==0){
      v_h_grs32_s3a[i]->GetXaxis()->SetTitle("energy [keV]");
      v_h_grs32_s3a[i]->GetYaxis()->SetTitle("counts / keV");
      v_h_grs32_s3a[i]->GetYaxis()->SetRangeUser(1., v_h_grs32_s3a[i]->GetMaximum()*10);
      v_h_grs32_s3a[i]->Draw("hist");
    }else v_h_grs32_s3a[i]->Draw("same hist");
  }
  leg_grs32_s3a->Draw("same");

  l_334->Draw("same");
  l_406->Draw("same");
  l_439->Draw("same");

  cc_grs32_s3a->SaveAs("fig/cc_grs32_s3a.png");

  //
  // ge ring4 sector1, s3 all
  TCanvas *cc_grs41_s3a = new TCanvas("cc_grs41_s3a","",900,500);
  cc_grs41_s3a->cd();
  cc_grs41_s3a->SetLogy();

  TLegend *leg_grs41_s3a = new TLegend(0.75, 0.38, 0.88, 0.92);
  for(int i=0;i<n_files;++i){
    leg_grs41_s3a->AddEntry(v_h_grs41_s3a[i], v_filename[i].substr(4,3).c_str());
    v_h_grs41_s3a[i]->SetTitle("ge ring4 sector1, different s3 beta pars");
    v_h_grs41_s3a[i]->SetLineColor(colors[i]);
    v_h_grs41_s3a[i]->SetLineWidth(2);
    v_h_grs41_s3a[i]->GetXaxis()->SetRangeUser(300,500);
    if(i==0){
      v_h_grs41_s3a[i]->GetXaxis()->SetTitle("energy [keV]");
      v_h_grs41_s3a[i]->GetYaxis()->SetTitle("counts / keV");
      v_h_grs41_s3a[i]->GetYaxis()->SetRangeUser(1., 100000);
      v_h_grs41_s3a[i]->Draw("hist");
    }else v_h_grs41_s3a[i]->Draw("same hist");
  }
  leg_grs41_s3a->Draw("same");

  l_334->Draw("same");
  l_406->Draw("same");
  l_439->Draw("same");

  cc_grs41_s3a->SaveAs("fig/cc_grs41_s3a.png");

  //
  // ge ring4 sector2, s3 all
  TCanvas *cc_grs42_s3a = new TCanvas("cc_grs42_s3a","",900,500);
  cc_grs42_s3a->cd();
  cc_grs42_s3a->SetLogy();

  TLegend *leg_grs42_s3a = new TLegend(0.75, 0.38, 0.88, 0.92);
  for(int i=0;i<n_files;++i){
    leg_grs42_s3a->AddEntry(v_h_grs42_s3a[i], v_filename[i].substr(4,3).c_str());
    v_h_grs42_s3a[i]->SetTitle("ge ring4 sector2, s3 differnet ring, different s3 beta pars");
    v_h_grs42_s3a[i]->SetLineColor(colors[i]);
    v_h_grs42_s3a[i]->SetLineWidth(2);
    v_h_grs42_s3a[i]->GetXaxis()->SetRangeUser(300,500);
    if(i==0){
      v_h_grs42_s3a[i]->GetXaxis()->SetTitle("energy [keV]");
      v_h_grs42_s3a[i]->GetYaxis()->SetTitle("counts / keV");
      v_h_grs42_s3a[i]->GetYaxis()->SetRangeUser(1., 100000);
      v_h_grs42_s3a[i]->Draw("hist");
    }else v_h_grs42_s3a[i]->Draw("same hist");
  }
  leg_grs42_s3a->Draw("same");

  l_334->Draw("same");
  l_406->Draw("same");
  l_439->Draw("same");

  cc_grs42_s3a->SaveAs("fig/cc_grs42_s3a.png");

  //
  // ge ring4 sector3, s3 all
  TCanvas *cc_grs43_s3a = new TCanvas("cc_grs43_s3a","",900,500);
  cc_grs43_s3a->cd();
  cc_grs43_s3a->SetLogy();

  TLegend *leg_grs43_s3a = new TLegend(0.75, 0.38, 0.88, 0.92);
  for(int i=0;i<n_files;++i){
    leg_grs43_s3a->AddEntry(v_h_grs43_s3a[i], v_filename[i].substr(4,3).c_str());
    v_h_grs43_s3a[i]->SetTitle("ge ring4 sector3, s3 differnet ring, different s3 beta pars");
    v_h_grs43_s3a[i]->SetLineColor(colors[i]);
    v_h_grs43_s3a[i]->SetLineWidth(2);
    v_h_grs43_s3a[i]->GetXaxis()->SetRangeUser(300,500);
    if(i==0){
      v_h_grs43_s3a[i]->GetXaxis()->SetTitle("energy [keV]");
      v_h_grs43_s3a[i]->GetYaxis()->SetTitle("counts / keV");
      v_h_grs43_s3a[i]->GetYaxis()->SetRangeUser(1., 100000);
      v_h_grs43_s3a[i]->Draw("hist");
    }else v_h_grs43_s3a[i]->Draw("same hist");
  }
  leg_grs43_s3a->Draw("same");

  l_334->Draw("same");
  l_406->Draw("same");
  l_439->Draw("same");

  cc_grs43_s3a->SaveAs("fig/cc_grs43_s3a.png");

  //
  // ge ring5 sector2, s3 all
  TCanvas *cc_grs52_s3a = new TCanvas("cc_grs52_s3a","",900,500);
  cc_grs52_s3a->cd();
  cc_grs52_s3a->SetLogy();

  TLegend *leg_grs52_s3a = new TLegend(0.75, 0.38, 0.88, 0.92);
  for(int i=0;i<n_files;++i){
    leg_grs52_s3a->AddEntry(v_h_grs52_s3a[i], v_filename[i].substr(4,3).c_str());
    v_h_grs52_s3a[i]->SetTitle("ge ring5 sector2, different s3 beta pars");
    v_h_grs52_s3a[i]->SetLineColor(colors[i]);
    v_h_grs52_s3a[i]->SetLineWidth(2);
    v_h_grs52_s3a[i]->GetXaxis()->SetRangeUser(300,500);
    if(i==0){
      v_h_grs52_s3a[i]->GetXaxis()->SetTitle("energy [keV]");
      v_h_grs52_s3a[i]->GetYaxis()->SetTitle("counts / keV");
      v_h_grs52_s3a[i]->GetYaxis()->SetRangeUser(1., 100000);
      v_h_grs52_s3a[i]->Draw("hist");
    }else v_h_grs52_s3a[i]->Draw("same hist");
  }
  leg_grs52_s3a->Draw("same");

  l_334->Draw("same");
  l_406->Draw("same");
  l_439->Draw("same");

  cc_grs52_s3a->SaveAs("fig/cc_grs52_s3a.png");

  //
  // ge ring5 sector3, s3 all
  TCanvas *cc_grs53_s3a = new TCanvas("cc_grs53_s3a","",900,500);
  cc_grs53_s3a->cd();
  cc_grs53_s3a->SetLogy();

  TLegend *leg_grs53_s3a = new TLegend(0.75, 0.38, 0.88, 0.92);
  for(int i=0;i<n_files;++i){
    leg_grs53_s3a->AddEntry(v_h_grs53_s3a[i], v_filename[i].substr(4,3).c_str());
    v_h_grs53_s3a[i]->SetTitle("ge ring5 sector3, different s3 beta pars");
    v_h_grs53_s3a[i]->SetLineColor(colors[i]);
    v_h_grs53_s3a[i]->SetLineWidth(2);
    v_h_grs53_s3a[i]->GetXaxis()->SetRangeUser(300,500);
    if(i==0){
      v_h_grs53_s3a[i]->GetXaxis()->SetTitle("energy [keV]");
      v_h_grs53_s3a[i]->GetYaxis()->SetTitle("counts / keV");
      v_h_grs53_s3a[i]->GetYaxis()->SetRangeUser(1., 100000);
      v_h_grs53_s3a[i]->Draw("hist");
    }else v_h_grs53_s3a[i]->Draw("same hist");
  }
  leg_grs53_s3a->Draw("same");

  l_334->Draw("same");
  l_406->Draw("same");
  l_439->Draw("same");

  cc_grs53_s3a->SaveAs("fig/cc_grs53_s3a.png");
}
