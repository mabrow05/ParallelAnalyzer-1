// Usage: root[0] .x plot_gain_bismuth.C("runNumber")

void plot_gain_LED(TString runNumber)
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Style options
  gROOT->SetStyle("Plain");
  //gStyle->SetOptStat(11);
  gStyle->SetOptStat(0);
  gStyle->SetStatFontSize(0.030);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFontSize(0.05);
  //gStyle->SetTitleX(0.17);
  //gStyle->SetTitleAlign(13);
  gStyle->SetTitleOffset(0.90, "x");
  gStyle->SetTitleOffset(1.25, "y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetNdivisions(510,"X");
  gStyle->SetNdivisions(510,"Y");
  gStyle->SetPadLeftMargin(0.15); // 0.13
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadBottomMargin(0.15); // 0.12
  gStyle->SetLabelSize(0.052,"X");
  gStyle->SetLabelSize(0.052,"Y");
  gStyle->SetTitleSize(0.052,"X");
  gStyle->SetTitleSize(0.052,"Y");
  gROOT->ForceStyle();
  gStyle->SetErrorX(0);

  // Open input file
  TString filenameIn;
  TString inDir = TString(getenv("GAIN_LED"));
  filenameIn  = inDir + "/gain_LED_";
  filenameIn += runNumber;
  filenameIn += ".root";
  cout << "Processing ... " << filenameIn << endl;
  TFile *filein = new TFile(filenameIn);

  // Output file
  TString filenameOut;
  filenameOut  = inDir + "/gain_LED_";
  filenameOut += runNumber;
  filenameOut += ".pdf";

  TString filenameOutFirst;
  filenameOutFirst  = filenameOut;
  filenameOutFirst += "(";

  TString filenameOutLast;
  filenameOutLast  = filenameOut;
  filenameOutLast += ")";

  // East
  c1 = new TCanvas("c1","c1");
  c1->Divide(2,2);

  c1_1->cd();
  c1_1->SetLogy(0);
  hisE0->SetXTitle("East PMT 1");
  hisE0->GetXaxis()->CenterTitle();
  hisE0->GetXaxis()->SetRangeUser(500.0,4000.0);
  hisE0->Draw();

  c1_2->cd();
  c1_2->SetLogy(0);
  hisE1->SetXTitle("East PMT 2");
  hisE1->GetXaxis()->CenterTitle();
  hisE1->GetXaxis()->SetRangeUser(500.0,4000.0);
  hisE1->Draw();

  c1_3->cd();
  c1_3->SetLogy(0);
  hisE2->SetXTitle("East PMT 3");
  hisE2->GetXaxis()->CenterTitle();
  hisE2->GetXaxis()->SetRangeUser(500.0,4000.0);
  hisE2->Draw();

  c1_4->cd();
  c1_4->SetLogy(0);
  hisE3->SetXTitle("East PMT 4");
  hisE3->GetXaxis()->CenterTitle();
  hisE3->GetXaxis()->SetRangeUser(500.0,4000.0);
  hisE3->Draw();

  c1->Print(filenameOutFirst);

  // West
  c2 = new TCanvas("c2","c2");
  c2->Divide(2,2);

  c2_1->cd();
  c2_1->SetLogy(0);
  hisW0->SetXTitle("West PMT 1");
  hisW0->GetXaxis()->CenterTitle();
  hisW0->GetXaxis()->SetRangeUser(500.0,4000.0);
  hisW0->Draw();

  c2_2->cd();
  c2_2->SetLogy(0);
  hisW1->SetXTitle("West PMT 2");
  hisW1->GetXaxis()->CenterTitle();
  hisW1->GetXaxis()->SetRangeUser(500.0,4000.0);
  hisW1->Draw();

  c2_3->cd();
  c2_3->SetLogy(0);
  hisW2->SetXTitle("West PMT 3");
  hisW2->GetXaxis()->CenterTitle();
  hisW2->GetXaxis()->SetRangeUser(500.0,4000.0);
  hisW2->Draw();

  c2_4->cd();
  c2_4->SetLogy(0);
  hisW3->SetXTitle("West PMT 4");
  hisW3->GetXaxis()->CenterTitle();
  hisW3->GetXaxis()->SetRangeUser(500.0,4000.0);
  hisW3->Draw();

  c2->Print(filenameOutLast);

}
