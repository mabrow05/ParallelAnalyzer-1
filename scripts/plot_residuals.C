#include "../include/sourcePeaks.h"

void plot_residuals(Int_t calLow, Int_t calHigh)
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Style options
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(11);
  //gStyle->SetOptStat(0);
  //gStyle->SetStatFontSize(0.020);
  gStyle->SetOptFit(1011); // 1111
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

  Int_t calPeriodLow = calLow;
  Int_t calPeriodHigh = calHigh;

  // Setup output file for error Envelope
  ofstream errEnv;
  Char_t tempfile[200];
  sprintf(tempfile,"%s/error_envelope/error_envelope_calPeriod_%i.dat",getenv("ANALYSIS_CODE"),calPeriodHigh);
  
  errEnv.open(tempfile);
  
  // Read  data file
  char temp[500];
  sprintf(temp, "%s/residuals/global_residuals_Erecon_runPeriods_%i-%i.dat",getenv("ANALYSIS_CODE"), calPeriodLow, calPeriodHigh);
  
  ifstream file(temp);

  const size_t N = 1000;
  TString source[N];
  Int_t run[N];
  Double_t res[N];
  Double_t resCe[N], resIn[N], resSn[N], resBi1[N], resBi2[N];

  Int_t i = 0;
  Int_t nCe = 0;
  Int_t nIn = 0;
  Int_t nSn = 0;
  Int_t nBi1 = 0;
  Int_t nBi2 = 0;
  cout << " \"bad\" runs: \n";
  while (!file.eof()) {
    file >> source[i] >> run[i] >> res[i];
    if (source[i] == "Ce") {
      resCe[nCe] = res[i];
      nCe++;
      if (sqrt(res[i]*res[i])>0.05*peakCe) cout << run[i] << " " << source[i] << " " << res[i] << endl; 
    }
    if (source[i] == "In") {
      resIn[nIn] = res[i];
      nIn++;
      if (sqrt(res[i]*res[i])>0.05*peakIn) cout << run[i] << " " << source[i] << " " << res[i] << endl; 
    }
    if (source[i] == "Sn") {
      resSn[nSn] = res[i];
      nSn++;
      if (sqrt(res[i]*res[i])>0.05*peakSn) cout << run[i] << " " << source[i] << " " << res[i] << endl; 
    }
    if (source[i] == "Bi1") {
      resBi1[nBi1] = res[i];
      nBi1++;
      if (sqrt(res[i]*res[i])>0.05*peakBiHigh) cout << run[i] << " " << source[i] << " " << res[i] << endl; 
    }
    if (source[i] == "Bi2") {
      resBi2[nBi2] = res[i];
      nBi2++;
      if (sqrt(res[i]*res[i])>0.05*peakBiLow) cout << run[i] << " " << source[i] << " " << res[i] << endl; 
    }
    if (file.fail()) break;
    i++;
  }
  cout << nCe << " " << nIn << " " << nSn << " " << nBi1 << " " << nBi2 << endl;

  
  // Histograms
  
  TH1F *hisCe = new TH1F("hisCe", "", 60, -30.0, 30.0);
  TH1F *hisIn = new TH1F("hisIn", "", 60, -30.0, 30.0);
  TH1F *hisSn = new TH1F("hisSn", "", 30, -30.0, 30.0);
  TH1F *hisBi1 = new TH1F("hisBi1", "", 30, -60.0, 60.0);
  TH1F *hisBi2 = new TH1F("hisBi2", "", 30, -60.0, 60.0);

  // Fill histograms
  for (int j=0; j<nCe; j++) {
    hisCe->Fill(resCe[j]);
  }
  for (int j=0; j<nIn; j++) {
    hisIn->Fill(resIn[j]);
  }
  for (int j=0; j<nSn; j++) {
    hisSn->Fill(resSn[j]);
  }
  for (int j=0; j<nBi1; j++) {
    hisBi1->Fill(resBi1[j]);
  }
  for (int j=0; j<nBi2; j++) {
    hisBi2->Fill(resBi2[j]);
  }
  // Calculate mean and standard deviation
  double meanCe = 0.;
  for (int j=0; j<nCe; j++) {
    meanCe += resCe[j];
  }
  meanCe = meanCe / (double) nCe;

  double sigmaCe = 0.;
  for (int j=0; j<nCe; j++) {
    sigmaCe += (resCe[j] - meanCe)*(resCe[j] - meanCe);
  }
  sigmaCe = sigmaCe / (double) nCe;
  sigmaCe = sqrt( sigmaCe );

  cout << "meanCe = " << meanCe << endl;
  cout << "     sigma = " << sigmaCe << endl;
  errEnv << "meanCe = " << meanCe << endl;
  errEnv << "sigma = " << sigmaCe << endl;

  double meanIn = 0.;
  for (int j=0; j<nIn; j++) {
    meanIn += resIn[j];
  }
  meanIn = meanIn / (double) nIn;

  double sigmaIn = 0.;
  for (int j=0; j<nIn; j++) {
    sigmaIn += (resIn[j] - meanIn)*(resIn[j] - meanIn);
  }
  sigmaIn = sigmaIn / (double) nIn;
  sigmaIn = sqrt( sigmaIn );

  cout << "meanIn = " << meanIn << endl;
  cout << "     sigma = " << sigmaIn << endl;
  errEnv << "meanIn = " << meanIn << endl;
  errEnv << "sigma = " << sigmaIn << endl;

  double meanSn = 0.;
  for (int j=0; j<nSn; j++) {
    meanSn += resSn[j];
  }
  meanSn = meanSn / (double) nSn;

  double sigmaSn = 0.;
  for (int j=0; j<nSn; j++) {
    sigmaSn += (resSn[j] - meanSn)*(resSn[j] - meanSn);
  }
  sigmaSn = sigmaSn / (double) nSn;
  sigmaSn = sqrt( sigmaSn );

  cout << "meanSn = " << meanSn << endl;
  cout << "     sigma = " << sigmaSn << endl;
  errEnv << "meanSn = " << meanSn << endl;
  errEnv << "sigma = " << sigmaSn << endl;

  double meanBi1 = 0.;
  for (int j=0; j<nBi1; j++) {
    meanBi1 += resBi1[j];
  }
  meanBi1 = meanBi1 / (double) nBi1;

  double sigmaBi1 = 0.;
  for (int j=0; j<nBi1; j++) {
    sigmaBi1 += (resBi1[j] - meanBi1)*(resBi1[j] - meanBi1);
  }
  sigmaBi1 = sigmaBi1 / (double) nBi1;
  sigmaBi1 = sqrt( sigmaBi1 );

  cout << "meanBi1 = " << meanBi1 << endl;
  cout << "     sigma = " << sigmaBi1 << endl;
  errEnv << "meanBi1 = " << meanBi1 << endl;
  errEnv << "sigma = " << sigmaBi1 << endl;

  double meanBi2 = 0.;
  for (int j=0; j<nBi2; j++) {
    meanBi2 += resBi2[j];
  }
  meanBi2 = meanBi2 / (double) nBi2;
  
  double sigmaBi2 = 0.;
  for (int j=0; j<nBi2; j++) {
    sigmaBi2 += (resBi2[j] - meanBi2)*(resBi2[j] - meanBi2);
  }
  sigmaBi2 = sigmaBi2 / (double) nBi2;
  sigmaBi2 = sqrt( sigmaBi2 );
  
  cout << "meanBi2 = " << meanBi2 << endl;
  cout << "     sigma = " << sigmaBi2 << endl;
  errEnv << "meanBi2 = " << meanBi2 << endl;
  errEnv << "sigma = " << sigmaBi2 << endl;

  
  
  cout << endl << "Results of Gaussian Fits:\n";
  errEnv << endl << "Results of Gaussian Fits:\n";

  // Ce 
  c1 = new TCanvas("c1", "c1");
  c1->SetLogy(0);

  hisCe->SetXTitle(" Ce E_{Q} Error [keV]");
  hisCe->GetXaxis()->CenterTitle();
  hisCe->SetLineColor(1);
  hisCe->Draw();
  hisCe->Fit("gaus", "", "", -8.0, 8.0);
  cout << "meanCe = " << gaus->GetParameter(1) << endl;
  cout << "     sigma = " << gaus->GetParameter(2)  << endl;
  errEnv << "meanCe = " << gaus->GetParameter(1) << endl;
  errEnv << "sigma = " << gaus->GetParameter(2) << endl;

  //In 
  c2 = new TCanvas("c2", "c2");
  c2->SetLogy(0);

  hisIn->SetXTitle(" In E_{Q} Error [keV]");
  hisIn->GetXaxis()->CenterTitle();
  hisIn->SetLineColor(1);
  hisIn->Draw();
  hisIn->Fit("gaus", "", "", -8.0, 8.0);
  cout << "meanIn = " << gaus->GetParameter(1) << endl;
  cout << "     sigma = " << gaus->GetParameter(2)  << endl;
  errEnv << "meanIn = " << gaus->GetParameter(1) << endl;
  errEnv << "sigma = " << gaus->GetParameter(2) << endl;

  // Sn 
  c3 = new TCanvas("c3", "c3");
  c3->SetLogy(0);

  hisSn->SetXTitle(" Sn E_{Q} Error [keV]");
  hisSn->GetXaxis()->CenterTitle();
  hisSn->SetLineColor(1);
  hisSn->Draw();
  hisSn->Fit("gaus", "", "", -25.0, 25.0);
  cout << "meanSn = " << gaus->GetParameter(1) << endl;
  cout << "     sigma = " << gaus->GetParameter(2)  << endl;
  errEnv << "meanSn = " << gaus->GetParameter(1) << endl;
  errEnv << "sigma = " << gaus->GetParameter(2) << endl;

  // Bi1 
  c4 = new TCanvas("c4", "c4");
  c4->SetLogy(0);

  hisBi1->SetXTitle(" Bi1 E_{Q} Error [keV]");
  hisBi1->GetXaxis()->CenterTitle();
  hisBi1->SetLineColor(1);
  hisBi1->Draw();
  hisBi1->Fit("gaus", "","", -46., 46.);
  cout << "meanBi1 = " << gaus->GetParameter(1) << endl;
  cout << "     sigma = " << gaus->GetParameter(2)  << endl;
  errEnv << "meanBi1 = " << gaus->GetParameter(1) << endl;
  errEnv << "sigma = " << gaus->GetParameter(2) << endl;

  // Bi2 
  c5 = new TCanvas("c5", "c5");
  c5->SetLogy(0);

  hisBi2->SetXTitle(" Bi2 E_{Q} Error [keV]");
  hisBi2->GetXaxis()->CenterTitle();
  hisBi2->SetLineColor(1);
  hisBi2->Draw();
  hisBi2->Fit("gaus", "","", -46., 46.);

  cout << "meanBi2 = " << gaus->GetParameter(1) << endl;
  cout << "     sigma = " << gaus->GetParameter(2)  << endl;
  errEnv << "meanBi2 = " << gaus->GetParameter(1) << endl;
  errEnv << "sigma = " << gaus->GetParameter(2) << endl;
  
  /*
  // Ce West
  c6 = new TCanvas("c6", "c6");
  c6->SetLogy(0);

  hisCeWest->SetXTitle("West Ce E_{Q} Error [keV]");
  hisCeWest->GetXaxis()->CenterTitle();
  hisCeWest->SetLineColor(1);
  hisCeWest->Draw();
  hisCeWest->Fit("gaus", "", "", -8., 8.);
  cout << "meanCeWest = " << gaus->GetParameter(1) << endl;
  cout << "     sigma = " << gaus->GetParameter(2)  << endl;
  errEnv << "meanCeWest = " << gaus->GetParameter(1) << endl;
  errEnv << "sigma = " << gaus->GetParameter(2) << endl;

  // In West
  c7 = new TCanvas("c7", "c7");
  c7->SetLogy(0);

  hisInWest->SetXTitle("West In E_{Q} Error [keV]");
  hisInWest->GetXaxis()->CenterTitle();
  hisInWest->SetLineColor(1);
  hisInWest->Draw();
  hisInWest->Fit("gaus", "", "", -8., 8.);
  cout << "meanInWest = " << gaus->GetParameter(1) << endl;
  cout << "     sigma = " << gaus->GetParameter(2)  << endl;
  errEnv << "meanInWest = " << gaus->GetParameter(1) << endl;
  errEnv << "sigma = " << gaus->GetParameter(2) << endl;


  // Sn West
  c8 = new TCanvas("c8", "c8");
  c8->SetLogy(0);

  hisSnWest->SetXTitle("West Sn E_{Q} Error [keV]");
  hisSnWest->GetXaxis()->CenterTitle();
  hisSnWest->SetLineColor(1);
  hisSnWest->Draw();
  hisSnWest->Fit("gaus", "", "", -25., 25.);
  cout << "meanCeWest = " << gaus->GetParameter(1) << endl;
  cout << "     sigma = " << gaus->GetParameter(2)  << endl;
  errEnv << "meanCeWest = " << gaus->GetParameter(1) << endl;
  errEnv << "sigma = " << gaus->GetParameter(2) << endl;

  // Bi1 West
  c9 = new TCanvas("c9", "c9");
  c9->SetLogy(0);

  hisBi1West->SetXTitle("West Bi1 E_{Q} Error [keV]");
  hisBi1West->GetXaxis()->CenterTitle();
  hisBi1West->SetLineColor(1);
  hisBi1West->Draw();
  hisBi1West->Fit("gaus", "", "", -46., 46.);
  cout << "meanBi1West = " << gaus->GetParameter(1) << endl;
  cout << "     sigma = " << gaus->GetParameter(2)  << endl;
  errEnv << "meanBi1West = " << gaus->GetParameter(1) << endl;
  errEnv << "sigma = " << gaus->GetParameter(2) << endl;

  // Bi2 West
  c10 = new TCanvas("c10", "c10");
  c10->SetLogy(0);

  hisBi2West->SetXTitle("West Bi2 E_{Q} Error [keV]");
  hisBi2West->GetXaxis()->CenterTitle();
  hisBi2West->SetLineColor(1);
  hisBi2West->Draw();
  hisBi2West->Fit("gaus", "", "", -46., 46.);
  if (postReplayPass4) {
    cout << "meanBi2West = " << gaus->GetParameter(1) << endl;
    cout << "     sigma = " << gaus->GetParameter(2)  << endl;
    errEnv << "meanBi2West = " << gaus->GetParameter(1) << endl;
    errEnv << "sigma = " << gaus->GetParameter(2) << endl;
  }

  // Ce 
  c11= new TCanvas("c11", "c11");
  c11->SetLogy(0);

  hisCe->SetXTitle("Total Ce E_{Q} Error [keV]");
  hisCe->GetXaxis()->CenterTitle();
  hisCe->SetLineColor(1);
  hisCe->Draw();
  hisCe->Fit("gaus", "", "", -8.0, 8.0);
  cout << "meanCe = " << gaus->GetParameter(1) << endl;
  cout << "     sigma = " << gaus->GetParameter(2)  << endl;
  errEnv << "meanCe = " << gaus->GetParameter(1) << endl;
  errEnv << "sigma = " << gaus->GetParameter(2) << endl;

  // In 
  c12= new TCanvas("c12", "c12");
  c12->SetLogy(0);

  hisIn->SetXTitle("Total In E_{Q} Error [keV]");
  hisIn->GetXaxis()->CenterTitle();
  hisIn->SetLineColor(1);
  hisIn->Draw();
  hisIn->Fit("gaus", "", "", -10.0, 10.0);
  cout << "meanIn = " << gaus->GetParameter(1) << endl;
  cout << "     sigma = " << gaus->GetParameter(2)  << endl;
  errEnv << "meanIn = " << gaus->GetParameter(1) << endl;
  errEnv << "sigma = " << gaus->GetParameter(2) << endl;

  // Sn 
  c13 = new TCanvas("c13", "c13");
  c13->SetLogy(0);

  hisSn->SetXTitle("Total Sn E_{Q} Error [keV]");
  hisSn->GetXaxis()->CenterTitle();
  hisSn->SetLineColor(1);
  hisSn->Draw();
  hisSn->Fit("gaus", "", "", -25.0, 25.0);
  cout << "meanSn = " << gaus->GetParameter(1) << endl;
  cout << "     sigma = " << gaus->GetParameter(2)  << endl;
  errEnv << "meanSn = " << gaus->GetParameter(1) << endl;
  errEnv << "sigma = " << gaus->GetParameter(2) << endl;

  // Bi1 
  c14 = new TCanvas("c14", "c14");
  c14->SetLogy(0);

  hisBi1->SetXTitle("Total Bi1 E_{Q} Error [keV]");
  hisBi1->GetXaxis()->CenterTitle();
  hisBi1->SetLineColor(1);
  hisBi1->Draw();
  hisBi1->Fit("gaus", "", "", -46., 46.);
  cout << "meanBi1 = " << gaus->GetParameter(1) << endl;
  cout << "     sigma = " << gaus->GetParameter(2)  << endl;
  errEnv << "meanBi1 = " << gaus->GetParameter(1) << endl;
  errEnv << "sigma = " << gaus->GetParameter(2) << endl;

  // Bi2 
  c15 = new TCanvas("c15", "c15");
  c15->SetLogy(0);

  hisBi2->SetXTitle("Total Bi2 E_{Q} Error [keV]");
  hisBi2->GetXaxis()->CenterTitle();
  hisBi2->SetLineColor(1);
  hisBi2->Draw();
  hisBi2->Fit("gaus", "", "", -46., 46.);
  if (postReplayPass4) {
    cout << "meanBi2 = " << gaus->GetParameter(1) << endl;
    cout << "     sigma = " << gaus->GetParameter(2)  << endl;
    errEnv << "meanBi2 = " << gaus->GetParameter(1) << endl;
    errEnv << "sigma = " << gaus->GetParameter(2) << endl;
  }*/

  errEnv.close();
}
