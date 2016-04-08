

{
  TFile *f = new TFile("test.root","READ");
  TGraphErrors *uk = (TGraphErrors*)f->Get("EreconALL_uk");
  TGraphErrors *sim = (TGraphErrors*)f->Get("EreconALL_sim");
  TCanvas *c1 = new TCanvas("c1");
  uk->SetMarkerColor(kRed);
  uk->SetLineColor(kRed);
  uk->SetFillStyle(3002);
  uk->SetFillColor(kRed);
  sim->SetMarkerColor(kBlue);
  uk->Draw("AB");
  sim->Draw("P");
}
