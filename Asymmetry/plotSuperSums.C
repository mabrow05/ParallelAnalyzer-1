

{
  int octetStart=0;
  int octetEnd=59;
  
  TString fileBase = "superSumPlots/SuperSum_octets_0-59";
  //TString::Format("%s.root",fileBase)
  TFile *f = new TFile("superSumPlots/SuperSum_octets_0-59.root","READ");

  TString pdfFile_base = TString::Format("%s.pdf",fileBase);
  TGraphErrors *uk = (TGraphErrors*)f->Get("EreconALL_uk");
  TGraphErrors *sim = (TGraphErrors*)f->Get("EreconALL_sim");
  TCanvas *c1 = new TCanvas("c1");
  
  uk->SetMarkerColor(kBlue);
  //uk->SetMarkerStyle
  uk->SetLineColor(kBlue);
  uk->SetFillStyle(3002);
  uk->SetFillColor(kBlue);
  sim->SetMarkerColor(kRed);
  //sim->SetMarkerSize(1);
  sim->SetMarkerStyle(20);
  uk->GetXaxis()->SetRangeUser(0., 800.);
  sim->GetXaxis()->SetRangeUser(0., 800.);
  uk->Draw("AB");
  sim->Draw("P");
}
