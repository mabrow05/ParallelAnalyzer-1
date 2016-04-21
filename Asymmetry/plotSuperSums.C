

{
  int octetStart=0;
  int octetEnd=59;
  
  TString fileBase = "superSumPlots/SuperSum_octets_0-59";
  //TString::Format("%s.root",fileBase)
  TFile *f = new TFile("superSumPlots/SuperSum_octets_0-59.root","READ");

  TString pdfFile_base = TString::Format("%s.pdf",fileBase);
  TGraphErrors *ukALL = (TGraphErrors*)f->Get("EreconALL_uk");
  TGraphErrors *simALL = (TGraphErrors*)f->Get("EreconALL_sim");
  TCanvas *c1 = new TCanvas("c1");
  
  ukALL->SetMarkerColor(kBlue);
  //ukALL->SetMarkerStyle
  ukALL->SetLineColor(kBlue);
  ukALL->SetFillStyle(3002);
  ukALL->SetFillColor(kBlue);
  simALL->SetMarkerColor(kRed);
  //sim->SetMarkerSize(1);
  simALL->SetMarkerStyle(20);
  ukALL->GetXaxis()->SetRangeUser(0., 800.);
  simALL->GetXaxis()->SetRangeUser(0., 800.);
  ukALL->Draw("AB");
  simALL->Draw("P");

  TGraphErrors *uk0 = (TGraphErrors*)f->Get("Erecon0_uk");
  TGraphErrors *sim0 = (TGraphErrors*)f->Get("Erecon0_sim");
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
  sim0->GetXaxis()->SetRangeUser(0., 800.);
  uk0->Draw("AB");
  sim0->Draw("P");

  TGraphErrors *uk1 = (TGraphErrors*)f->Get("Erecon1_uk");
  TGraphErrors *sim1 = (TGraphErrors*)f->Get("Erecon1_sim");
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
  sim1->GetXaxis()->SetRangeUser(0., 800.);
  uk1->Draw("AB");
  sim1->Draw("P");

  TGraphErrors *uk23 = (TGraphErrors*)f->Get("Erecon23_uk");
  TGraphErrors *sim23 = (TGraphErrors*)f->Get("Erecon23_sim");
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
  sim23->GetXaxis()->SetRangeUser(0., 800.);
  uk23->Draw("AB");
  sim23->Draw("P");
}
