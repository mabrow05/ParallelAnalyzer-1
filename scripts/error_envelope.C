{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Style options
  gROOT->SetStyle("Plain");
  //gStyle->SetOptStat(11);
  gStyle->SetOptStat(0);
  gStyle->SetStatFontSize(0.020);
  gStyle->SetOptFit(0); // 1111
  gStyle->SetOptTitle(1);
  gStyle->SetTitleFontSize(0.05);
  //gStyle->SetTitleX(0.17);
  //gStyle->SetTitleAlign(13);
  gStyle->SetTitleOffset(1.20, "x");
  gStyle->SetTitleOffset(1.10, "y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetNdivisions(510,"X");
  gStyle->SetNdivisions(510,"Y");
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.12);

  // East residuals
  double EQ[3];
  EQ[0] = 98.2;
  EQ[1] = 331.2;
  EQ[2] = 928.0;
  double dEQ[3] = {};

  double resEast[3];
  resEast[0] = -0.2002;
  resEast[1] =  0.5212;
  resEast[2] = -.3176;

  double sigEast[3];
  sigEast[0] = 1.8833;
  sigEast[1] = 5.29;
  sigEast[2] = 10.157;

  // Plot
  c1 = new TCanvas("c1", "canvas");
  c1->SetLogy(0);

  TGraphErrors *gr1 = new TGraphErrors(3,EQ,resEast,dEQ,sigEast);
  gr1->SetTitle("");
  gr1->SetMarkerColor(1);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1.0);

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr1,"P");
  mg->Draw("A");
  mg->GetXaxis()->SetTitle("E_{Q} [keV]");
  mg->GetXaxis()->SetTitleOffset(1.2);
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->SetTitle("Calibration Residual [keV]");
  mg->GetYaxis()->SetTitleOffset(1.2);
  mg->GetYaxis()->CenterTitle();

  mg->GetXaxis()->SetLimits(0.0,1000.0);
  mg->SetMinimum(-30.0);
  mg->SetMaximum( 30.0);

  Int_t n = 2;
  Double_t x[n] = {0, 1000};
  Double_t y[n] = {0.0, 0.0};

  TGraph *gr0 = new TGraph(n,x,y);
  gr0->Draw("Same");
  gr0->SetLineWidth(2);
  gr0->SetLineColor(1);
  gr0->SetLineStyle(2);

}
