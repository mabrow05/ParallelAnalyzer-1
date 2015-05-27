

void MB_errorEnvelope(Int_t calLow, Int_t calHigh, Int_t pmt)
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

  Int_t PMT = pmt; //0->Average over PMTs; 1,2,3,4 -> single PMT

  Int_t calPeriodLow = calLow;
  Int_t calPeriodHigh = calHigh;

  // Setup output file for error Envelope
  ofstream errEnv;
  Char_t tempfile[200];
  if (calPeriodLow!=calPeriodHigh && !PMT) {
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
  if (calPeriodLow!=calPeriodHigh && !PMT) sprintf(tempEast, "../residuals/residuals_global_East_periods_%i-%i.dat", calPeriodLow, calPeriodHigh);
  else if (calPeriodLow==calPeriodHigh && !PMT) sprintf(tempEast, "../residuals/residuals_East_runPeriod_%i.dat", calPeriodLow);
  else if (calPeriodLow!=calPeriodHigh && PMT) sprintf(tempEast,"../residuals/residuals_global_East_periods_%i-%i_PMTE%i.dat", calPeriodLow, calPeriodHigh, PMT);
  else if (calPeriodLow==calPeriodHigh && PMT) sprintf(tempEast,"../residuals/residuals_East_runPeriod_PMTE%i.dat", calPeriodLow,PMT);
  ifstream fileEast(tempEast);

  const size_t N = 500;
  TString sourceEast[N];
  Int_t runEast[N];
  Double_t resEast[N];
  Double_t resCeEast[N], resSnEast[N], resBiEast[N];

  // E_Q peaks 
  Double_t peakCe = 98.2;
  Double_t peakSn = 331.2;
  Double_t peakBiLow = 443.0;
  Double_t peakBiHigh = 928.0;

  Int_t i = 0;
  Int_t nCeEast = 0;
  Int_t nSnEast = 0;
  Int_t nBiEast = 0;
  cout << "East \"bad\" runs: \n";
  while (!fileEast.eof()) {
    fileEast >> sourceEast[i] >> runEast[i] >> resEast[i];
    if (sourceEast[i] == "Ce_East") {
      resCeEast[nCeEast] = resEast[i];
      nCeEast++;
      if (sqrt(resEast[i]*resEast[i])>0.05*peakCe) cout << runEast[i] << " " << sourceEast[i] << " " << resEast[i] << endl; 
    }
    if (sourceEast[i] == "Sn_East") {
      resSnEast[nSnEast] = resEast[i];
      nSnEast++;
      if (sqrt(resEast[i]*resEast[i])>0.05*peakSn) cout << runEast[i] << " " << sourceEast[i] << " " << resEast[i] << endl; 
    }
    if (sourceEast[i] == "Bi_East") {
      resBiEast[nBiEast] = resEast[i];
      nBiEast++;
      if (sqrt(resEast[i]*resEast[i])>0.05*peakBiHigh) cout << runEast[i] << " " << sourceEast[i] << " " << resEast[i] << endl; 
    }
    if (fileEast.fail()) break;
    i++;
  }
  cout << nCeEast << " " << nSnEast << " " << nBiEast << endl;

  // Read West data file
  char tempWest[500];
  if (calPeriodLow!=calPeriodHigh && !PMT) sprintf(tempWest, "../residuals/residuals_global_West_periods_%i-%i.dat", calPeriodLow, calPeriodHigh);
  else if (calPeriodLow==calPeriodHigh && !PMT) sprintf(tempWest, "../residuals/residuals_West_runPeriod_%i.dat", calPeriodLow);
  else if (calPeriodLow!=calPeriodHigh && PMT) sprintf(tempWest,"../residuals/residuals_global_West_periods_%i-%i_PMTW%i.dat", calPeriodLow, calPeriodHigh, PMT);
  else if (calPeriodLow==calPeriodHigh && PMT) sprintf(tempWest,"../residuals/residuals_West_runPeriod_PMTW%i.dat", calPeriodLow,PMT);
  ifstream fileWest(tempWest);

  TString sourceWest[N];
  Int_t runWest[N];
  Double_t resWest[N];
  Double_t resCeWest[N], resSnWest[N], resBiWest[N];

  Int_t i = 0;
  Int_t nCeWest = 0;
  Int_t nSnWest = 0;
  Int_t nBiWest = 0;
  
  cout << "West \"bad\" runs: \n";
  while (!fileWest.eof()) {
    fileWest >> sourceWest[i] >> runWest[i] >> resWest[i];
    if (sourceWest[i] == "Ce_West") {
      resCeWest[nCeWest] = resWest[i];
      nCeWest++;
      if (sqrt(resWest[i]*resWest[i])>0.05*peakCe) cout << runWest[i] << " " << sourceWest[i] << " " << resWest[i] << endl; 
    }
    if (sourceWest[i] == "Sn_West") {
      resSnWest[nSnWest] = resWest[i];
      nSnWest++;
      if (sqrt(resWest[i]*resWest[i])>0.05*peakSn) cout << runWest[i] << " " << sourceWest[i] << " " << resWest[i] << endl; 
    }
    if (sourceWest[i] == "Bi_West") {
      resBiWest[nBiWest] = resWest[i];
      nBiWest++;
      if (sqrt(resWest[i]*resWest[i])>0.05*peakBiHigh) cout << runWest[i] << " " << sourceWest[i] << " " << resWest[i] << endl; 
    }
    if (fileWest.fail()) break;
    i++;
  }
  cout << nCeWest << " " << nSnWest << " " << nBiWest << endl;

  // Histograms
  TH1F *hisCeEast = new TH1F("hisCeEast", "", 60, -30.0, 30.0);
  TH1F *hisSnEast = new TH1F("hisSnEast", "", 30, -30.0, 30.0);
  TH1F *hisBiEast = new TH1F("hisBiEast", "", 30, -60.0, 60.0);
  TH1F *hisCeWest = new TH1F("hisCeWest", "", 60, -30.0, 30.0);
  TH1F *hisSnWest = new TH1F("hisSnWest", "", 30, -30.0, 30.0);
  TH1F *hisBiWest = new TH1F("hisBiWest", "", 30, -60.0, 60.0);

  // Fill histograms
  for (int j=0; j<nCeEast; j++) {
    hisCeEast->Fill(resCeEast[j]);
  }
  for (int j=0; j<nSnEast; j++) {
    hisSnEast->Fill(resSnEast[j]);
  }
  for (int j=0; j<nBiEast; j++) {
    hisBiEast->Fill(resBiEast[j]);
  }
  for (int j=0; j<nCeWest; j++) {
    hisCeWest->Fill(resCeWest[j]);
  }
  for (int j=0; j<nSnWest; j++) {
    hisSnWest->Fill(resSnWest[j]);
  }
  for (int j=0; j<nBiWest; j++) {
    hisBiWest->Fill(resBiWest[j]);
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

  double meanBiEast = 0.;
  for (int j=0; j<nBiEast; j++) {
    meanBiEast += resBiEast[j];
  }
  meanBiEast = meanBiEast / (double) nBiEast;

  double sigmaBiEast = 0.;
  for (int j=0; j<nBiEast; j++) {
    sigmaBiEast += (resBiEast[j] - meanBiEast)*(resBiEast[j] - meanBiEast);
  }
  sigmaBiEast = sigmaBiEast / (double) nBiEast;
  sigmaBiEast = sqrt( sigmaBiEast );

  cout << "meanBiEast = " << meanBiEast << endl;
  cout << "     sigma = " << sigmaBiEast << endl;
  errEnv << "meanBiEast = " << meanBiEast << endl;
  errEnv << "sigma = " << sigmaBiEast << endl;

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

  double meanBiWest = 0.;
  for (int j=0; j<nBiWest; j++) {
    meanBiWest += resBiWest[j];
  }
  meanBiWest = meanBiWest / (double) nBiWest;

  double sigmaBiWest = 0.;
  for (int j=0; j<nBiWest; j++) {
    sigmaBiWest += (resBiWest[j] - meanBiWest)*(resBiWest[j] - meanBiWest);
  }
  sigmaBiWest = sigmaBiWest / (double) nBiWest;
  sigmaBiWest = sqrt( sigmaBiWest );

  cout << "meanBiWest = " << meanBiWest << endl;
  cout << "     sigma = " << sigmaBiWest << endl;
  errEnv << "meanBiWest = " << meanBiWest << endl;
  errEnv << "sigma = " << sigmaBiWest << endl;

  errEnv.close();

  // Ce East
  c1 = new TCanvas("c1", "c1");
  c1->SetLogy(0);

  hisCeEast->SetXTitle("East Ce E_{Q} Error [keV]");
  hisCeEast->GetXaxis()->CenterTitle();
  hisCeEast->SetLineColor(1);
  hisCeEast->Draw();
  //hisCeEast->Fit("gaus", "", "", -5.0, 5.0);

  // Sn East
  c2 = new TCanvas("c2", "c2");
  c2->SetLogy(0);

  hisSnEast->SetXTitle("East Sn E_{Q} Error [keV]");
  hisSnEast->GetXaxis()->CenterTitle();
  hisSnEast->SetLineColor(1);
  hisSnEast->Draw();
  //hisSnEast->Fit("gaus", "", "", -10.0, 10.0);

  // Bi East
  c3 = new TCanvas("c3", "c3");
  c3->SetLogy(0);

  hisBiEast->SetXTitle("East Bi E_{Q} Error [keV]");
  hisBiEast->GetXaxis()->CenterTitle();
  hisBiEast->SetLineColor(1);
  hisBiEast->Draw();
  //hisBiEast->Fit("gaus", "");

  // Ce West
  c4 = new TCanvas("c4", "c4");
  c4->SetLogy(0);

  hisCeWest->SetXTitle("West Ce E_{Q} Error [keV]");
  hisCeWest->GetXaxis()->CenterTitle();
  hisCeWest->SetLineColor(1);
  hisCeWest->Draw();
  //hisCeWest->Fit("gaus", "");

  // Sn West
  c5 = new TCanvas("c5", "c5");
  c5->SetLogy(0);

  hisSnWest->SetXTitle("West Sn E_{Q} Error [keV]");
  hisSnWest->GetXaxis()->CenterTitle();
  hisSnWest->SetLineColor(1);
  hisSnWest->Draw();
  //hisSnWest->Fit("gaus", "");

  // Bi West
  c6 = new TCanvas("c6", "c6");
  c6->SetLogy(0);

  hisBiWest->SetXTitle("West Bi E_{Q} Error [keV]");
  hisBiWest->GetXaxis()->CenterTitle();
  hisBiWest->SetLineColor(1);
  hisBiWest->Draw();
  //hisBiWest->Fit("gaus", "");

}
