

{
  gStyle->SetOptStat(0);

  int octetStart=0;
  int octetEnd=59;
  
  //TString fileBase = "superSumPlots/SuperSum_octets_0-59";
  //TString::Format("%s.root",fileBase)
  TFile *f = new TFile("superSumPlots/SuperSum_octets_0-59_norm_60.0-450.0.root","READ");

  //TString pdfFile_base = TString::Format("%s.pdf",fileBase);
  TH1D *ukALL = (TH1D*)f->Get("EreconALL_uk");
  TH1D *simALL = (TH1D*)f->Get("EreconALL_sim");

  //Residuals...
  TCanvas *r = new TCanvas("r");
  Int_t nBins = ukALL->GetNbinsX();
  Double_t min = ukALL->GetXaxis()->GetBinLowEdge(ukALL->GetXaxis()->GetFirst());
  Double_t max = ukALL->GetXaxis()->GetBinUpEdge(ukALL->GetXaxis()->GetLast());
  TH1D *resid = new TH1D("resid","Residuals: MC-Data", nBins, min, max);
  resid->Add(simALL,ukALL,100,-100);
  resid->GetXaxis()->SetRangeUser(0., 800.);
  resid->GetYaxis()->SetTitle("event rate (mHz/keV)");
  resid->Draw();
  r->Update();
  TLine *l = new TLine(min, 0., 800., 0.);
  l->SetLineStyle(8);
  l->Draw();

  TCanvas *c1 = new TCanvas("c1");
  
  ukALL->SetMarkerColor(kBlue);
  ukALL->SetMarkerStyle(22)
  ukALL->SetLineColor(kBlue);
  //ukALL->SetFillStyle(3002);
  //ukALL->SetFillColor(kBlue);
  //ukALL->SetLineWidth(3);
  simALL->SetMarkerColor(kRed);
  simALL->SetLineColor(kRed);
  //simALL->SetLineWidth(3);
  //sim->SetMarkerSize(1);
  simALL->SetMarkerStyle(20);
  ukALL->GetXaxis()->SetRangeUser(0., 800.);
  ukALL->GetYaxis()->SetTitle("event rate (mHz/keV)");
  ukALL->Scale(100.);
  simALL->Scale(100.);
  simALL->GetXaxis()->SetRangeUser(0., 800.);
  ukALL->Draw("");
  simALL->Draw("SAME");
  
  TH1D *uk0 = (TH1D*)f->Get("Erecon0_uk");
  TH1D *sim0 = (TH1D*)f->Get("Erecon0_sim");
  TCanvas *c2 = new TCanvas("c2");
  
  uk0->SetMarkerColor(kBlue);
  //uk0->SetMarkerStyle
  uk0->SetLineColor(kBlue);
  uk0->SetFillStyle(3002);
  uk0->SetFillColor(kBlue);
  sim0->SetMarkerColor(kRed);
  //sim->SetMarkerSize(1);
  sim0->SetMarkerStyle(20);
  uk0->GetXaxis()->SetRangeUser(0., 800.);
  uk0->GetYaxis()->SetTitle("event rate (mHz/keV)");
  uk0->Scale(100.);
  sim0->Scale(100.);
  sim0->GetXaxis()->SetRangeUser(0., 800.);
  uk0->Draw();
  sim0->Draw("SAME");

  TH1D *uk1 = (TH1D*)f->Get("Erecon1_uk");
  TH1D *sim1 = (TH1D*)f->Get("Erecon1_sim");
  TCanvas *c3 = new TCanvas("c3");
  
  uk1->SetMarkerColor(kBlue);
  //uk1->SetMarkerStyle
  uk1->SetLineColor(kBlue);
  uk1->SetFillStyle(3002);
  uk1->SetFillColor(kBlue);
  sim1->SetMarkerColor(kRed);
  //sim->SetMarkerSize(1);
  sim1->SetMarkerStyle(20);
  uk1->GetXaxis()->SetRangeUser(0., 800.);
  uk1->GetYaxis()->SetTitle("event rate (mHz/keV)");
  uk1->Scale(100.);
  sim1->Scale(100.);
  sim1->GetXaxis()->SetRangeUser(0., 800.);
  uk1->Draw();
  sim1->Draw("SAME");

  TH1D *uk23 = (TH1D*)f->Get("Erecon23_uk");
  TH1D *sim23 = (TH1D*)f->Get("Erecon23_sim");
  TCanvas *c4 = new TCanvas("c4");
  
  uk23->SetMarkerColor(kBlue);
  //uk23->SetMarkerStyle
  uk23->SetLineColor(kBlue);
  uk23->SetFillStyle(3002);
  uk23->SetFillColor(kBlue);
  sim23->SetMarkerColor(kRed);
  //sim->SetMarkerSize(1);
  sim23->SetMarkerStyle(20);
  uk23->GetXaxis()->SetRangeUser(0., 800.);
  uk23->GetYaxis()->SetTitle("event rate (mHz/keV)");
  uk23->Scale(100.);
  sim23->Scale(100.);
  sim23->GetXaxis()->SetRangeUser(0., 800.);
  uk23->Draw();
  sim23->Draw("SAME");

  
  /*for (int i=1; i<=nBins; i++) {
    Double_t BinContent = simALL->GetBinContent(i)-ukALL->GetBinContent(i);
    resid->SetBinContent(i,BinContent);
    Double_t BinError = sqrt(simALL->GetBinError(i)*simALL->GetBinError(i)+ukALL->GetBinError(i)*ukALL->GetBinError(i);
    resid->SetBinError(i,sqrt(simALL->GetBinError(i)*/
}
