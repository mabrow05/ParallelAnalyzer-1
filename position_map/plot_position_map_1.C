{

  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Style options
  gROOT->SetStyle("Plain");
  Int_t palette[10];
  palette[0] = 0;
  palette[1] = 19;
  palette[2] = 18;
  palette[3] = 17;
  palette[4] = 16;
  palette[5] = 15;
  palette[6] = 14;
  palette[7] = 13;
  palette[8] = 12;
  palette[9] = 1;
  gStyle->SetPalette(10,palette);  // z-axis color scale for 2D histograms
  //gStyle->SetPalette(2);  // z-axis color scale for 2D histograms
  //gStyle->SetOptStat(11);
  gStyle->SetOptStat(0);
  gStyle->SetStatFontSize(0.030);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFontSize(0.05);
  //gStyle->SetTitleX(0.17);
  //gStyle->SetTitleAlign(13);
  gStyle->SetTitleOffset(0.80, "x");
  gStyle->SetTitleOffset(1.30, "y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetNdivisions(510,"X");
  gStyle->SetNdivisions(510,"Y");
  gStyle->SetNdivisions(9,"Z");
  gStyle->SetPadLeftMargin(0.13); // 0.13
  gStyle->SetPadRightMargin(0.15); // 0.04
  gStyle->SetPadBottomMargin(0.37); // 0.30
  gStyle->SetLabelSize(0.045, "X");
  gStyle->SetLabelSize(0.045, "Y");
  gStyle->SetLabelSize(0.045, "Z");
  gStyle->SetLabelOffset(0.00, "X");
  gStyle->SetLabelOffset(0.01, "Y");
  gStyle->SetTitleSize(0.050, "X");
  gStyle->SetTitleSize(0.050, "Y");
  gROOT->ForceStyle();

  // Geometry
  int nPosBinsX = 11;
  int nPosBinsY = 11;
  double xBinWidth = 10.0;
  double yBinWidth = 10.0;
  double xBinLower[nPosBinsX];
  double xBinUpper[nPosBinsX];
  double xBinCenter[nPosBinsX];
  double yBinLower[nPosBinsY];
  double yBinUpper[nPosBinsY];
  double yBinCenter[nPosBinsY];

  for (int k=0; k<nPosBinsX; k++) {
    xBinLower[k]     = -55.0 + ((double) k)*xBinWidth;
    xBinUpper[k]     = -55.0 + ((double) k)*xBinWidth + xBinWidth;
    xBinCenter[k]    = (xBinLower[k] + xBinUpper[k])/2.;
  }
  for (int k=0; k<nPosBinsY; k++) {
    yBinLower[k]     = -55.0 + ((double) k)*yBinWidth;
    yBinUpper[k]     = -55.0 + ((double) k)*yBinWidth + yBinWidth;
    yBinCenter[k]    = (yBinLower[k] + yBinUpper[k])/2.;
  }

  // Histograms
  TH2F *hisE0 = new TH2F("E0", "", nPosBinsX,-55.0,55.0, nPosBinsY,-55.0,55.0);
  TH2F *hisE1 = new TH2F("E1", "", nPosBinsX,-55.0,55.0, nPosBinsY,-55.0,55.0);
  TH2F *hisE2 = new TH2F("E2", "", nPosBinsX,-55.0,55.0, nPosBinsY,-55.0,55.0);
  TH2F *hisE3 = new TH2F("E3", "", nPosBinsX,-55.0,55.0, nPosBinsY,-55.0,55.0);
  TH2F *hisW0 = new TH2F("W0", "", nPosBinsX,-55.0,55.0, nPosBinsY,-55.0,55.0);
  TH2F *hisW1 = new TH2F("W1", "", nPosBinsX,-55.0,55.0, nPosBinsY,-55.0,55.0);
  TH2F *hisW2 = new TH2F("W2", "", nPosBinsX,-55.0,55.0, nPosBinsY,-55.0,55.0);
  TH2F *hisW3 = new TH2F("W3", "", nPosBinsX,-55.0,55.0, nPosBinsY,-55.0,55.0);

  // Read position map
  int nPMT = 8;
  char tempIn[500];
  sprintf(tempIn, "/extern/plaster/ucna/Analyzer/position_map/position_map_1.dat");
  ifstream fileIn(tempIn);

  double posMap[nPMT][nPosBinsX][nPosBinsY];

  for (int i=0; i<nPosBinsX; i++) {
    for (int j=0; j<nPosBinsY; j++) {
      fileIn >> posMap[0][i][j]
             >> posMap[1][i][j]
             >> posMap[2][i][j]
             >> posMap[3][i][j]
             >> posMap[4][i][j]
             >> posMap[5][i][j]
             >> posMap[6][i][j]
	     >> posMap[7][i][j];
    }
  }

  // Fill histograms
  for (int i=0; i<nPosBinsX; i++) {
    for (int j=0; j<nPosBinsY; j++) {
      hisE0->Fill(xBinCenter[i],yBinCenter[j],posMap[0][i][j]);
      hisE1->Fill(xBinCenter[i],yBinCenter[j],posMap[1][i][j]);
      hisE2->Fill(xBinCenter[i],yBinCenter[j],posMap[2][i][j]);
      hisE3->Fill(xBinCenter[i],yBinCenter[j],posMap[3][i][j]);
      hisW0->Fill(xBinCenter[i],yBinCenter[j],posMap[4][i][j]);
      hisW1->Fill(xBinCenter[i],yBinCenter[j],posMap[5][i][j]);
      hisW2->Fill(xBinCenter[i],yBinCenter[j],posMap[6][i][j]);
      hisW3->Fill(xBinCenter[i],yBinCenter[j],posMap[7][i][j]);
    }
  }

  // Output PDF file
  TString filenameOut;
  filenameOut  = "position_map_1.pdf";

  TString filenameOutFirst;
  filenameOutFirst  = filenameOut;
  filenameOutFirst += "(";

  TString filenameOutLast;
  filenameOutLast  = filenameOut;
  filenameOutLast += ")";

  c0 = new TCanvas("c0", "canvas");
  c0->Divide(2,1);
  c0->SetLogy(0);

  // East 0
  c0_1->cd();
  hisE0->Draw("colz");
  hisE0->SetXTitle("x [mm]");
  hisE0->SetYTitle("y [mm]");
  hisE0->GetXaxis()->CenterTitle();
  hisE0->GetYaxis()->CenterTitle();
  hisE0->GetYaxis()->CenterTitle();
  hisE0.SetAxisRange(0.5,1.5,"Z");

  /*
  TEllipse *ell = new TEllipse(0,0,45,45);
  ell->SetFillStyle(0);
  ell->SetLineColor(1);
  ell->SetLineStyle(2);
  ell->SetLineWidth(3);
  ell->Draw("same");
  */

  Double_t x1_text =  40;
  Double_t y1_text =  45;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.052);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("E0");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  // East 1
  c0_2->cd();
  hisE1->Draw("colz");
  hisE1->SetXTitle("x [mm]");
  hisE1->SetYTitle("y [mm]");
  hisE1->GetXaxis()->CenterTitle();
  hisE1->GetYaxis()->CenterTitle();
  hisE1->GetYaxis()->CenterTitle();
  hisE1.SetAxisRange(0.5,1.5,"Z");

  /*
  TEllipse *ell = new TEllipse(0,0,45,45);
  ell->SetFillStyle(0);
  ell->SetLineColor(1);
  ell->SetLineStyle(2);
  ell->SetLineWidth(3);
  ell->Draw("same");
  */

  Double_t x1_text = -45;
  Double_t y1_text =  45;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.052);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("E1");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  c0->Print(filenameOutFirst);

  c1 = new TCanvas("c1", "canvas");
  c1->Divide(2,1);
  c1->SetLogy(0);

  // East 2
  c1_1->cd();
  hisE2->Draw("colz");
  hisE2->SetXTitle("x [mm]");
  hisE2->SetYTitle("y [mm]");
  hisE2->GetXaxis()->CenterTitle();
  hisE2->GetYaxis()->CenterTitle();
  hisE2->GetYaxis()->CenterTitle();
  hisE2.SetAxisRange(0.5,1.5,"Z");

  /*
  TEllipse *ell = new TEllipse(0,0,45,45);
  ell->SetFillStyle(0);
  ell->SetLineColor(1);
  ell->SetLineStyle(2);
  ell->SetLineWidth(3);
  ell->Draw("same");
  */

  Double_t x1_text = -45;
  Double_t y1_text = -45;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.052);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("E2");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  // East 3
  c1_2->cd();
  hisE3->Draw("colz");
  hisE3->SetXTitle("x [mm]");
  hisE3->SetYTitle("y [mm]");
  hisE3->GetXaxis()->CenterTitle();
  hisE3->GetYaxis()->CenterTitle();
  hisE3->GetYaxis()->CenterTitle();
  hisE3.SetAxisRange(0.5,1.5,"Z");

  /*
  TEllipse *ell = new TEllipse(0,0,45,45);
  ell->SetFillStyle(0);
  ell->SetLineColor(1);
  ell->SetLineStyle(2);
  ell->SetLineWidth(3);
  ell->Draw("same");
  */

  Double_t x1_text =  40;
  Double_t y1_text = -45;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.052);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("E3");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  c1->Print(filenameOut);

  c2 = new TCanvas("c2", "canvas");
  c2->Divide(2,1);
  c2->SetLogy(0);

  // West 0
  c2_1->cd();
  hisW0->Draw("colz");
  hisW0->SetXTitle("x [mm]");
  hisW0->SetYTitle("y [mm]");
  hisW0->GetXaxis()->CenterTitle();
  hisW0->GetYaxis()->CenterTitle();
  hisW0->GetYaxis()->CenterTitle();
  hisW0.SetAxisRange(0.5,1.5,"Z");

  /*
  TEllipse *ell = new TEllipse(0,0,45,45);
  ell->SetFillStyle(0);
  ell->SetLineColor(1);
  ell->SetLineStyle(2);
  ell->SetLineWidth(3);
  ell->Draw("same");
  */

  Double_t x1_text = -45;
  Double_t y1_text = -45;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.052);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("W0");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  // West 1
  c2_2->cd();
  hisW1->Draw("colz");
  hisW1->SetXTitle("x [mm]");
  hisW1->SetYTitle("y [mm]");
  hisW1->GetXaxis()->CenterTitle();
  hisW1->GetYaxis()->CenterTitle();
  hisW1->GetYaxis()->CenterTitle();
  hisW1.SetAxisRange(0.5,1.5,"Z");

  /*
  TEllipse *ell = new TEllipse(0,0,45,45);
  ell->SetFillStyle(0);
  ell->SetLineColor(1);
  ell->SetLineStyle(2);
  ell->SetLineWidth(3);
  ell->Draw("same");
  */

  Double_t x1_text =  40;
  Double_t y1_text = -45;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.052);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("W1");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  c2->Print(filenameOut);

  c3 = new TCanvas("c3", "canvas");
  c3->Divide(2,1);
  c3->SetLogy(0);

  // West 2
  c3_1->cd();
  hisW2->Draw("colz");
  hisW2->SetXTitle("x [mm]");
  hisW2->SetYTitle("y [mm]");
  hisW2->GetXaxis()->CenterTitle();
  hisW2->GetYaxis()->CenterTitle();
  hisW2->GetYaxis()->CenterTitle();
  hisW2.SetAxisRange(0.5,1.5,"Z");

  /*
  TEllipse *ell = new TEllipse(0,0,45,45);
  ell->SetFillStyle(0);
  ell->SetLineColor(1);
  ell->SetLineStyle(2);
  ell->SetLineWidth(3);
  ell->Draw("same");
  */

  Double_t x1_text =  40;
  Double_t y1_text =  45;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.052);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("W2");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  // West 3
  c3_2->cd();
  hisW3->Draw("colz");
  hisW3->SetXTitle("x [mm]");
  hisW3->SetYTitle("y [mm]");
  hisW3->GetXaxis()->CenterTitle();
  hisW3->GetYaxis()->CenterTitle();
  hisW3->GetYaxis()->CenterTitle();
  hisW3.SetAxisRange(0.5,1.5,"Z");

  /*
  TEllipse *ell = new TEllipse(0,0,45,45);
  ell->SetFillStyle(0);
  ell->SetLineColor(1);
  ell->SetLineStyle(2);
  ell->SetLineWidth(3);
  ell->Draw("same");
  */

  Double_t x1_text = -45;
  Double_t y1_text =  45;

  TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
  pt1->SetTextSize(0.052);
  pt1->SetTextColor(1);
  pt1->SetTextAlign(12);
  pt1->AddText("W3");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  pt1->SetFillColor(0);
  pt1->Draw();

  c3->Print(filenameOutLast);

  // Redraw axis covered up by gray band
  gPad->RedrawAxis();
}
