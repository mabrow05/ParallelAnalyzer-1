#include "../include/sourcePeaks.h"

void MB_errorEnvelope(Int_t calLow, Int_t calHigh, Int_t pmt, bool postReplayPass4)
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

  Int_t PMT = pmt; //0->Average over PMTs; 1,2,3,4 -> single PMT

  Int_t calPeriodLow = calLow;
  Int_t calPeriodHigh = calHigh;

  // Setup output file for error Envelope
  ofstream errEnv;
  Char_t tempfile[200];
  if (postReplayPass4) {
    sprintf(tempfile,"../error_envelope/error_envelope_PostReplayPass4_PMTbyPMT_calPeriods_%i-%i.dat",calPeriodLow,calPeriodHigh);}
  else if (calPeriodLow!=calPeriodHigh && !PMT) {
    sprintf(tempfile,"../error_envelope/error_envelope_calPeriods_%i-%i.dat",calPeriodLow,calPeriodHigh);}
  else if (calPeriodLow!=calPeriodHigh && PMT) {
    sprintf(tempfile,"../error_envelope/error_envelope_calPeriods_%i-%i_PMT%i.dat",calPeriodLow,calPeriodHigh,PMT);}
  else if (calPeriodLow==calPeriodHigh && !PMT) {
    sprintf(tempfile,"../error_envelope/error_envelope_calPeriod_%i.dat",calPeriodHigh);}
  else if (calPeriodLow==calPeriodHigh && PMT) {
    sprintf(tempfile,"../error_envelope/error_envelope_calPeriod_%i_PMT%i.dat",calPeriodLow,PMT);}

  errEnv.open(tempfile);

  // Read East data file
  char tempEast[500];
  if (postReplayPass4) sprintf(tempEast, "../residuals/residuals_global_EvisPMTbyPMT_East_periods_%i-%i.dat", calPeriodLow, calPeriodHigh);
  else if (calPeriodLow!=calPeriodHigh && !PMT) sprintf(tempEast, "../residuals/residuals_global_East_periods_%i-%i.dat", calPeriodLow, calPeriodHigh);
  else if (calPeriodLow==calPeriodHigh && !PMT) sprintf(tempEast, "../residuals/residuals_East_runPeriod_%i.dat", calPeriodLow);
  else if (calPeriodLow!=calPeriodHigh && PMT) sprintf(tempEast,"../residuals/residuals_global_East_periods_%i-%i_PMTE%i.dat", calPeriodLow, calPeriodHigh, PMT);
  else if (calPeriodLow==calPeriodHigh && PMT) sprintf(tempEast,"../residuals/residuals_East_runPeriod_%i_PMTE%i.dat", calPeriodLow,PMT);
  ifstream fileEast(tempEast);

  const size_t N = 1000;
  TString sourceEast[N];
  Int_t runEast[N];
  Double_t resEast[N];
  Double_t resCeEast[N], resSnEast[N], resBi1East[N], resBi2East[N];

  Int_t i = 0;
  Int_t nCeEast = 0;
  Int_t nSnEast = 0;
  Int_t nBi1East = 0;
  Int_t nBi2East = 0;
  cout << "East \"bad\" runs: \n";
  while (fileEast >> sourceEast[i] >> runEast[i] >> resEast[i]) {
    if (sourceEast[i] == "Ce_East") {
      resCeEast[nCeEast] = resEast[i];
      nCeEast++;
      if (sqrt(resEast[i]*resEast[i])>0.05*peakCe_EQ) cout << runEast[i] << " " << sourceEast[i] << " " << resEast[i] << endl; 
    }
    if (sourceEast[i] == "Sn_East") {
      resSnEast[nSnEast] = resEast[i];
      nSnEast++;
      if (sqrt(resEast[i]*resEast[i])>0.05*peakSn_EQ) cout << runEast[i] << " " << sourceEast[i] << " " << resEast[i] << endl; 
    }
    if (sourceEast[i] == "Bi1_East") {
      resBi1East[nBi1East] = resEast[i];
      nBi1East++;
      if (sqrt(resEast[i]*resEast[i])>0.05*peakBiHigh_EQ) cout << runEast[i] << " " << sourceEast[i] << " " << resEast[i] << endl; 
    }
    if (sourceEast[i] == "Bi2_East") {
      resBi2East[nBi2East] = resEast[i];
      nBi2East++;
      if (sqrt(resEast[i]*resEast[i])>0.05*peakBiLow_EQ) cout << runEast[i] << " " << sourceEast[i] << " " << resEast[i] << endl; 
    }
    if (fileEast.fail()) break;
    i++;
  }
  cout << nCeEast << " " << nSnEast << " " << nBi1East << " " << nBi2East << endl;

  // Read West data file
  char tempWest[500];
  if (postReplayPass4) sprintf(tempWest, "../residuals/residuals_global_EvisPMTbyPMT_West_periods_%i-%i.dat", calPeriodLow, calPeriodHigh);
  else if (calPeriodLow!=calPeriodHigh && !PMT) sprintf(tempWest, "../residuals/residuals_global_West_periods_%i-%i.dat", calPeriodLow, calPeriodHigh);
  else if (calPeriodLow==calPeriodHigh && !PMT) sprintf(tempWest, "../residuals/residuals_West_runPeriod_%i.dat", calPeriodLow);
  else if (calPeriodLow!=calPeriodHigh && PMT) sprintf(tempWest,"../residuals/residuals_global_West_periods_%i-%i_PMTW%i.dat", calPeriodLow, calPeriodHigh, PMT);
  else if (calPeriodLow==calPeriodHigh && PMT) sprintf(tempWest,"../residuals/residuals_West_runPeriod_%i_PMTW%i.dat", calPeriodLow,PMT);
  ifstream fileWest(tempWest);

  TString sourceWest[N];
  Int_t runWest[N];
  Double_t resWest[N];
  Double_t resCeWest[N], resSnWest[N], resBi1West[N], resBi2West[N];

  Int_t i = 0;
  Int_t nCeWest = 0;
  Int_t nSnWest = 0;
  Int_t nBi1West = 0;
  Int_t nBi2West = 0;
  
  cout << "West \"bad\" runs: \n";
  while (fileWest >> sourceWest[i] >> runWest[i] >> resWest[i]) {
    if (sourceWest[i] == "Ce_West") {
      resCeWest[nCeWest] = resWest[i];
      nCeWest++;
      if (sqrt(resWest[i]*resWest[i])>0.05*peakCe_EQ) cout << runWest[i] << " " << sourceWest[i] << " " << resWest[i] << endl; 
    }
    if (sourceWest[i] == "Sn_West") {
      resSnWest[nSnWest] = resWest[i];
      nSnWest++;
      if (sqrt(resWest[i]*resWest[i])>0.05*peakSn_EQ) cout << runWest[i] << " " << sourceWest[i] << " " << resWest[i] << endl; 
    }    
    if (sourceWest[i] == "Bi1_West") {
      resBi1West[nBi1West] = resWest[i];
      nBi1West++;
      if (sqrt(resWest[i]*resWest[i])>0.05*peakBiHigh_EQ) cout << runWest[i] << " " << sourceWest[i] << " " << resWest[i] << endl; 
    }
    if (sourceWest[i] == "Bi2_West") {
      resBi2West[nBi2West] = resWest[i];
      nBi2West++;
      if (sqrt(resWest[i]*resWest[i])>0.05*peakBiLow_EQ) cout << runWest[i] << " " << sourceWest[i] << " " << resWest[i] << endl; 
    }
    if (fileWest.fail()) break;
    i++;
  }
  cout << nCeWest << " " << nSnWest << " " << nBi1West << " " << nBi2West << endl;

  // Histograms
  TH1F *hisCeEast = new TH1F("hisCeEast", "", 60, -30.0, 30.0);
  TH1F *hisSnEast = new TH1F("hisSnEast", "", 30, -30.0, 30.0);
  TH1F *hisBi1East = new TH1F("hisBi1East", "", 30, -60.0, 60.0);
  TH1F *hisBi2East = new TH1F("hisBi2East", "", 30, -60.0, 60.0);
  TH1F *hisCeWest = new TH1F("hisCeWest", "", 60, -30.0, 30.0);
  TH1F *hisSnWest = new TH1F("hisSnWest", "", 30, -30.0, 30.0);
  TH1F *hisBi1West = new TH1F("hisBi1West", "", 30, -60.0, 60.0);
  TH1F *hisBi2West = new TH1F("hisBi2West", "", 30, -60.0, 60.0);
  TH1F *hisCe = new TH1F("hisCe", "", 60, -30.0, 30.0);
  TH1F *hisSn = new TH1F("hisSn", "", 30, -30.0, 30.0);
  TH1F *hisBi1 = new TH1F("hisBi1", "", 30, -60.0, 60.0);
  TH1F *hisBi2 = new TH1F("hisBi2", "", 30, -60.0, 60.0);

  // Fill histograms
  for (int j=0; j<nCeEast; j++) {
    hisCeEast->Fill(resCeEast[j]);
    hisCe->Fill(resCeEast[j]);
  }
  for (int j=0; j<nSnEast; j++) {
    hisSnEast->Fill(resSnEast[j]);
    hisSn->Fill(resSnEast[j]);
  }
  for (int j=0; j<nBi1East; j++) {
    hisBi1East->Fill(resBi1East[j]);
    hisBi1->Fill(resBi1East[j]);
  }
  for (int j=0; j<nBi2East; j++) {
    hisBi2East->Fill(resBi2East[j]);
    hisBi2->Fill(resBi2East[j]);
  }
  for (int j=0; j<nCeWest; j++) {
    hisCeWest->Fill(resCeWest[j]);
    hisCe->Fill(resCeWest[j]);
  }
  for (int j=0; j<nSnWest; j++) {
    hisSnWest->Fill(resSnWest[j]);
    hisSn->Fill(resSnWest[j]);
  }
  for (int j=0; j<nBi1West; j++) {
    hisBi1West->Fill(resBi1West[j]);
    hisBi1->Fill(resBi1West[j]);
  }
 for (int j=0; j<nBi2West; j++) {
    hisBi2West->Fill(resBi2West[j]);
    hisBi2->Fill(resBi2West[j]);
  }
  // Calculate mean and standard deviation
  double meanCeEast = 0.;
  for (int j=0; j<nCeEast; j++) {
    meanCeEast += resCeEast[j];
  }
  meanCeEast = meanCeEast / (double) nCeEast;

  double sigmaCeEast = 0.;
  for (int j=0; j<nCeEast; j++) {
    sigmaCeEast += (resCeEast[j] - meanCeEast)*(resCeEast[j] - meanCeEast);
  }
  sigmaCeEast = sigmaCeEast / (double) nCeEast;
  sigmaCeEast = sqrt( sigmaCeEast );

  cout << "meanCeEast = " << meanCeEast << endl;
  cout << "     sigma = " << sigmaCeEast << endl;
  errEnv << "meanCeEast = " << meanCeEast << endl;
  errEnv << "sigma = " << sigmaCeEast << endl;

  double meanSnEast = 0.;
  for (int j=0; j<nSnEast; j++) {
    meanSnEast += resSnEast[j];
  }
  meanSnEast = meanSnEast / (double) nSnEast;

  double sigmaSnEast = 0.;
  for (int j=0; j<nSnEast; j++) {
    sigmaSnEast += (resSnEast[j] - meanSnEast)*(resSnEast[j] - meanSnEast);
  }
  sigmaSnEast = sigmaSnEast / (double) nSnEast;
  sigmaSnEast = sqrt( sigmaSnEast );

  cout << "meanSnEast = " << meanSnEast << endl;
  cout << "     sigma = " << sigmaSnEast << endl;
  errEnv << "meanSnEast = " << meanSnEast << endl;
  errEnv << "sigma = " << sigmaSnEast << endl;

  double meanBi1East = 0.;
  for (int j=0; j<nBi1East; j++) {
    meanBi1East += resBi1East[j];
  }
  meanBi1East = meanBi1East / (double) nBi1East;

  double sigmaBi1East = 0.;
  for (int j=0; j<nBi1East; j++) {
    sigmaBi1East += (resBi1East[j] - meanBi1East)*(resBi1East[j] - meanBi1East);
  }
  sigmaBi1East = sigmaBi1East / (double) nBi1East;
  sigmaBi1East = sqrt( sigmaBi1East );

  cout << "meanBi1East = " << meanBi1East << endl;
  cout << "     sigma = " << sigmaBi1East << endl;
  errEnv << "meanBi1East = " << meanBi1East << endl;
  errEnv << "sigma = " << sigmaBi1East << endl;

  if (postReplayPass4) {
    double meanBi2East = 0.;
    for (int j=0; j<nBi2East; j++) {
      meanBi2East += resBi2East[j];
    }
    meanBi2East = meanBi2East / (double) nBi2East;
    
    double sigmaBi2East = 0.;
    for (int j=0; j<nBi2East; j++) {
      sigmaBi2East += (resBi2East[j] - meanBi2East)*(resBi2East[j] - meanBi2East);
    }
    sigmaBi2East = sigmaBi2East / (double) nBi2East;
    sigmaBi2East = sqrt( sigmaBi2East );
    
    cout << "meanBi2East = " << meanBi2East << endl;
    cout << "     sigma = " << sigmaBi2East << endl;
    errEnv << "meanBi2East = " << meanBi2East << endl;
    errEnv << "sigma = " << sigmaBi2East << endl;
  }
  // Calculate mean and standard deviation
  double meanCeWest = 0.;
  for (int j=0; j<nCeWest; j++) {
    meanCeWest += resCeWest[j];
  }
  meanCeWest = meanCeWest / (double) nCeWest;

  double sigmaCeWest = 0.;
  for (int j=0; j<nCeWest; j++) {
    sigmaCeWest += (resCeWest[j] - meanCeWest)*(resCeWest[j] - meanCeWest);
  }
  sigmaCeWest = sigmaCeWest / (double) nCeWest;
  sigmaCeWest = sqrt( sigmaCeWest );

  cout << "meanCeWest = " << meanCeWest << endl;
  cout << "     sigma = " << sigmaCeWest << endl;
  errEnv << "meanCeWest = " << meanCeWest << endl;
  errEnv << "sigma = " << sigmaCeWest << endl;

  double meanSnWest = 0.;
  for (int j=0; j<nSnWest; j++) {
    meanSnWest += resSnWest[j];
  }
  meanSnWest = meanSnWest / (double) nSnWest;

  double sigmaSnWest = 0.;
  for (int j=0; j<nSnWest; j++) {
    sigmaSnWest += (resSnWest[j] - meanSnWest)*(resSnWest[j] - meanSnWest);
  }
  sigmaSnWest = sigmaSnWest / (double) nSnWest;
  sigmaSnWest = sqrt( sigmaSnWest );

  cout << "meanSnWest = " << meanSnWest << endl;
  cout << "     sigma = " << sigmaSnWest << endl;
  errEnv << "meanSnWest = " << meanSnWest << endl;
  errEnv << "sigma = " << sigmaSnWest << endl;

  double meanBi1West = 0.;
  for (int j=0; j<nBi1West; j++) {
    meanBi1West += resBi1West[j];
  }
  meanBi1West = meanBi1West / (double) nBi1West;

  double sigmaBi1West = 0.;
  for (int j=0; j<nBi1West; j++) {
    sigmaBi1West += (resBi1West[j] - meanBi1West)*(resBi1West[j] - meanBi1West);
  }
  sigmaBi1West = sigmaBi1West / (double) nBi1West;
  sigmaBi1West = sqrt( sigmaBi1West );

  cout << "meanBi1West = " << meanBi1West << endl;
  cout << "     sigma = " << sigmaBi1West << endl;
  errEnv << "meanBi1West = " << meanBi1West << endl;
  errEnv << "sigma = " << sigmaBi1West << endl;

  if (postReplayPass4) {
    double meanBi2West = 0.;
    for (int j=0; j<nBi2West; j++) {
      meanBi2West += resBi2West[j];
    }
    meanBi2West = meanBi2West / (double) nBi2West;
    
    double sigmaBi2West = 0.;
    for (int j=0; j<nBi2West; j++) {
      sigmaBi2West += (resBi2West[j] - meanBi2West)*(resBi2West[j] - meanBi2West);
    }
    sigmaBi2West = sigmaBi2West / (double) nBi2West;
    sigmaBi2West = sqrt( sigmaBi2West );
    
    cout << "meanBi2West = " << meanBi2West << endl;
    cout << "     sigma = " << sigmaBi2West << endl;
    errEnv << "meanBi2West = " << meanBi2West << endl;
    errEnv << "sigma = " << sigmaBi2West << endl;
  }
  //Calculate total mean and rms
  double meanCe = 0.;
  for (int j=0; j<nCeEast; j++) {
    meanCe += resCeEast[j];
  }
  for (int j=0; j<nCeWest; j++) {
    meanCe += resCeWest[j];
  }
  meanCe = meanCe / (double) (nCeEast+nCeWest);

  double sigmaCe = 0.;
  for (int j=0; j<nCeEast; j++) {
    sigmaCe += (resCeEast[j] - meanCe)*(resCeEast[j] - meanCe);
  }
  for (int j=0; j<nCeWest; j++) {
    sigmaCe += (resCeWest[j] - meanCe)*(resCeWest[j] - meanCe);
  }
  sigmaCe = sigmaCe / (double) (nCeEast+nCeWest);
  sigmaCe = sqrt( sigmaCe );

  cout << "meanCe = " << meanCe << endl;
  cout << "     sigma = " << sigmaCe << endl;
  errEnv << "meanCe = " << meanCe << endl;
  errEnv << "sigma = " << sigmaCe << endl;


  double meanSn = 0.;
  for (int j=0; j<nSnEast; j++) {
    meanSn += resSnEast[j];
  }
  for (int j=0; j<nSnWest; j++) {
    meanSn += resSnWest[j];
  }
  meanSn = meanSn / (double) (nSnEast+nSnWest);

  double sigmaSn = 0.;
  for (int j=0; j<nSnEast; j++) {
    sigmaSn += (resSnEast[j] - meanSn)*(resSnEast[j] - meanSn);
  }
  for (int j=0; j<nSnWest; j++) {
    sigmaSn += (resSnWest[j] - meanSn)*(resSnWest[j] - meanSn);
  }
  sigmaSn = sigmaSn / (double) (nSnEast+nSnWest);
  sigmaSn = sqrt( sigmaSn );

  cout << "meanSn = " << meanSn << endl;
  cout << "     sigma = " << sigmaSn << endl;
  errEnv << "meanSn = " << meanSn << endl;
  errEnv << "sigma = " << sigmaSn << endl;


  double meanBi1 = 0.;
  for (int j=0; j<nBi1East; j++) {
    meanBi1 += resBi1East[j];
  }
  for (int j=0; j<nBi1West; j++) {
    meanBi1 += resBi1West[j];
  }
  meanBi1 = meanBi1 / (double) (nBi1East+nBi1West);

  double sigmaBi1 = 0.;
  for (int j=0; j<nBi1East; j++) {
    sigmaBi1 += (resBi1East[j] - meanBi1)*(resBi1East[j] - meanBi1);
  }
  for (int j=0; j<nBi1West; j++) {
    sigmaBi1 += (resBi1West[j] - meanBi1)*(resBi1West[j] - meanBi1);
  }
  sigmaBi1 = sigmaBi1 / (double) (nBi1East+nBi1West);
  sigmaBi1 = sqrt( sigmaBi1 );

  cout << "meanBi1 = " << meanBi1 << endl;
  cout << "     sigma = " << sigmaBi1 << endl;
  errEnv << "meanBi1 = " << meanBi1 << endl;
  errEnv << "sigma = " << sigmaBi1 << endl;

  if (postReplayPass4) {
    double meanBi2 = 0.;
    for (int j=0; j<nBi2East; j++) {
      meanBi2 += resBi2East[j];
    }
    for (int j=0; j<nBi2West; j++) {
      meanBi2 += resBi2West[j];
    }
    meanBi2 = meanBi2 / (double) (nBi2East+nBi2West);
    
    double sigmaBi2 = 0.;
    for (int j=0; j<nBi2East; j++) {
      sigmaBi2 += (resBi2East[j] - meanBi2)*(resBi2East[j] - meanBi2);
    }
    for (int j=0; j<nBi2West; j++) {
      sigmaBi2 += (resBi2West[j] - meanBi2)*(resBi2West[j] - meanBi2);
    }
    sigmaBi2 = sigmaBi2 / (double) (nBi2East+nBi2West);
    sigmaBi2 = sqrt( sigmaBi2 );
    
    cout << "meanBi2 = " << meanBi2 << endl;
    cout << "     sigma = " << sigmaBi2 << endl;
    errEnv << "meanBi2 = " << meanBi2 << endl;
    errEnv << "sigma = " << sigmaBi2 << endl;
  }

  cout << endl << "Results of Gaussian Fits:\n";
  errEnv << endl << "Results of Gaussian Fits:\n";

  // Ce East
  c1 = new TCanvas("c1", "c1");
  c1->SetLogy(0);

  hisCeEast->SetXTitle("East Ce E_{Q} Error [keV]");
  hisCeEast->GetXaxis()->CenterTitle();
  hisCeEast->SetLineColor(1);
  hisCeEast->Draw();
  hisCeEast->Fit("gaus", "", "", -8.0, 8.0);
  cout << "meanCeEast = " << gaus->GetParameter(1) << endl;
  cout << "     sigma = " << gaus->GetParameter(2)  << endl;
  errEnv << "meanCeEast = " << gaus->GetParameter(1) << endl;
  errEnv << "sigma = " << gaus->GetParameter(2) << endl;

  // Sn East
  c2 = new TCanvas("c2", "c2");
  c2->SetLogy(0);

  hisSnEast->SetXTitle("East Sn E_{Q} Error [keV]");
  hisSnEast->GetXaxis()->CenterTitle();
  hisSnEast->SetLineColor(1);
  hisSnEast->Draw();
  hisSnEast->Fit("gaus", "", "", -25.0, 25.0);
  cout << "meanSnEast = " << gaus->GetParameter(1) << endl;
  cout << "     sigma = " << gaus->GetParameter(2)  << endl;
  errEnv << "meanSnEast = " << gaus->GetParameter(1) << endl;
  errEnv << "sigma = " << gaus->GetParameter(2) << endl;

  // Bi1 East
  c3 = new TCanvas("c3", "c3");
  c3->SetLogy(0);

  hisBi1East->SetXTitle("East Bi1 E_{Q} Error [keV]");
  hisBi1East->GetXaxis()->CenterTitle();
  hisBi1East->SetLineColor(1);
  hisBi1East->Draw();
  hisBi1East->Fit("gaus", "","", -46., 46.);
  cout << "meanBi1East = " << gaus->GetParameter(1) << endl;
  cout << "     sigma = " << gaus->GetParameter(2)  << endl;
  errEnv << "meanBi1East = " << gaus->GetParameter(1) << endl;
  errEnv << "sigma = " << gaus->GetParameter(2) << endl;

  // Bi2 East
  c4 = new TCanvas("c4", "c4");
  c4->SetLogy(0);

  hisBi2East->SetXTitle("East Bi2 E_{Q} Error [keV]");
  hisBi2East->GetXaxis()->CenterTitle();
  hisBi2East->SetLineColor(1);
  hisBi2East->Draw();
  hisBi2East->Fit("gaus", "","", -46., 46.);
  if (postReplayPass4) {
    cout << "meanBi2East = " << gaus->GetParameter(1) << endl;
    cout << "     sigma = " << gaus->GetParameter(2)  << endl;
    errEnv << "meanBi2East = " << gaus->GetParameter(1) << endl;
    errEnv << "sigma = " << gaus->GetParameter(2) << endl;
  }

  // Ce West
  c5 = new TCanvas("c5", "c5");
  c5->SetLogy(0);

  hisCeWest->SetXTitle("West Ce E_{Q} Error [keV]");
  hisCeWest->GetXaxis()->CenterTitle();
  hisCeWest->SetLineColor(1);
  hisCeWest->Draw();
  hisCeWest->Fit("gaus", "", "", -8., 8.);
  cout << "meanCeWest = " << gaus->GetParameter(1) << endl;
  cout << "     sigma = " << gaus->GetParameter(2)  << endl;
  errEnv << "meanCeWest = " << gaus->GetParameter(1) << endl;
  errEnv << "sigma = " << gaus->GetParameter(2) << endl;

  // Sn West
  c6 = new TCanvas("c6", "c6");
  c6->SetLogy(0);

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
  c7 = new TCanvas("c7", "c7");
  c7->SetLogy(0);

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
  c8 = new TCanvas("c8", "c8");
  c8->SetLogy(0);

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
  c9 = new TCanvas("c9", "c9");
  c9->SetLogy(0);

  hisCe->SetXTitle("Total Ce E_{Q} Error [keV]");
  hisCe->GetXaxis()->CenterTitle();
  hisCe->SetLineColor(1);
  hisCe->Draw();
  hisCe->Fit("gaus", "", "", -8.0, 8.0);
  cout << "meanCe = " << gaus->GetParameter(1) << endl;
  cout << "     sigma = " << gaus->GetParameter(2)  << endl;
  errEnv << "meanCe = " << gaus->GetParameter(1) << endl;
  errEnv << "sigma = " << gaus->GetParameter(2) << endl;

  // Sn 
  c10 = new TCanvas("c10", "c10");
  c10->SetLogy(0);

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
  c11 = new TCanvas("c11", "c11");
  c11->SetLogy(0);

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
  c12 = new TCanvas("c12", "c12");
  c12->SetLogy(0);

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
  }

  errEnv.close();
}
