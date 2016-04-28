

{
  gStyle->SetOptStat(0);

  int octetStart=0;
  int octetEnd=59;
  
  TString fileBase = "superSumPlots/SuperSum_octets_0-59_norm_0.0-1200.0";
  TString fileName = TString::Format("%s.root",fileBase.Data());
  TFile *f = new TFile(fileName,"READ");

  
  TH1D *ukALL = (TH1D*)f->Get("EreconALL_uk");
  TH1D *simALL = (TH1D*)f->Get("EreconALL_sim");

  //Residuals...
  TCanvas *c1 = new TCanvas("c1","c1", 1600., 600.);
  c1->Divide(2,1);
  c1->cd(2);
  Int_t nBins = ukALL->GetNbinsX();
  Double_t min = ukALL->GetXaxis()->GetBinLowEdge(ukALL->GetXaxis()->GetFirst());
  Double_t max = ukALL->GetXaxis()->GetBinUpEdge(ukALL->GetXaxis()->GetLast());
  TH1D *resid = new TH1D("resid","Residuals: MC-Data", nBins, min, max);
  resid->Add(simALL,ukALL,100,-100);
  resid->GetXaxis()->SetRangeUser(0., 1200.);
  resid->GetYaxis()->SetTitle("event rate (mHz/keV)");
  resid->SetLineWidth(2);
  resid->Draw();
  c1->Update();
  TLine *l = new TLine(min, 0., 1200., 0.);
  l->SetLineStyle(8);
  l->Draw();

  c1->cd(1);
  
  ukALL->SetMarkerColor(kBlue);
  ukALL->SetLineColor(kBlue);
  ukALL->SetFillStyle(3002);
  ukALL->SetFillColor(kBlue);
  ukALL->SetLineWidth(3);
  simALL->SetMarkerColor(kRed);
  simALL->SetLineColor(kRed);
  simALL->SetMarkerSize(0.75);
  simALL->SetMarkerStyle(20);
  ukALL->GetXaxis()->SetRangeUser(0., 1200.);
  ukALL->GetYaxis()->SetTitle("event rate (mHz/keV)");
  ukALL->Scale(100.);
  simALL->Scale(100.);
  simALL->GetXaxis()->SetRangeUser(0., 1200.);
  ukALL->Draw("BARE0");
  simALL->Draw("SAMEE0");


  TCanvas *c2 = new TCanvas("c2", "c2", 1600., 600.);
  c2->Divide(3,1);
  c2->cd(1);
  
  TH1D *uk0 = (TH1D*)f->Get("Erecon0_uk");
  TH1D *sim0 = (TH1D*)f->Get("Erecon0_sim");

  uk0->SetMarkerColor(kBlue);
  uk0->SetLineColor(kBlue);
  uk0->SetFillStyle(3002);
  uk0->SetFillColor(kBlue);
  uk0->SetLineWidth(3);
  sim0->SetMarkerColor(kRed);
  sim0->SetLineColor(kRed);
  sim0->SetMarkerSize(0.75);
  sim0->SetMarkerStyle(20);
  uk0->GetXaxis()->SetRangeUser(0., 1200.);
  uk0->GetYaxis()->SetTitle("event rate (mHz/keV)");
  uk0->Scale(100.);
  sim0->Scale(100.);
  sim0->GetXaxis()->SetRangeUser(0., 1200.);
  uk0->Draw("BARE0");
  sim0->Draw("SAMEE0");
  
  c2->cd(2);
  TH1D *uk1 = (TH1D*)f->Get("Erecon1_uk");
  TH1D *sim1 = (TH1D*)f->Get("Erecon1_sim");

  uk1->SetMarkerColor(kBlue);
  uk1->SetLineColor(kBlue);
  uk1->SetFillStyle(3002);
  uk1->SetFillColor(kBlue);
  uk1->SetLineWidth(3);
  sim1->SetMarkerColor(kRed);
  sim1->SetLineColor(kRed);
  sim1->SetMarkerSize(0.75);
  sim1->SetMarkerStyle(20);
  uk1->GetXaxis()->SetRangeUser(0., 1200.);
  uk1->GetYaxis()->SetTitle("event rate (mHz/keV)");
  uk1->Scale(100.);
  sim1->Scale(100.);
  sim1->GetXaxis()->SetRangeUser(0., 1200.);
  uk1->Draw("BARE0");
  sim1->Draw("SAMEE0");
  

  c2->cd(3);

  TH1D *uk23 = (TH1D*)f->Get("Erecon23_uk");
  TH1D *sim23 = (TH1D*)f->Get("Erecon23_sim");
  
  uk23->SetMarkerColor(kBlue);
  uk23->SetLineColor(kBlue);
  uk23->SetFillStyle(3002);
  uk23->SetFillColor(kBlue);
  uk23->SetLineWidth(3);
  sim23->SetMarkerColor(kRed);
  sim23->SetLineColor(kRed);
  sim23->SetMarkerSize(0.75);
  sim23->SetMarkerStyle(20);
  uk23->GetXaxis()->SetRangeUser(0., 1200.);
  uk23->GetYaxis()->SetTitle("event rate (mHz/keV)");
  uk23->Scale(100.);
  sim23->Scale(100.);
  sim23->GetXaxis()->SetRangeUser(0., 1200.);
  uk23->Draw("BARE0");
  sim23->Draw("SAMEE0");


  TString pdfFile = TString::Format("%s.pdf",fileBase.Data());
  c1->Print(pdfFile);
  c1->Print(pdfFile);
  
}
