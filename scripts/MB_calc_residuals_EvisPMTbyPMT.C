#include <vector>
#include <algorithm>
#include <map>
#include "../include/sourcePeaks.h"

unsigned int find_vec_location_int(std::vector<Int_t> vec, Int_t val)
{
  UInt_t size = vec.size();
  //for (int i=0;i<size;i++){cout << vec[i] << endl;}

  std::vector<Int_t>::iterator it = std::find(vec.begin(),vec.end(),val);
  if (it!=vec.end()) {
    return it - vec.begin();
  }
  else cout << "Can't locate run " << val << " in PMT Quality list" << endl; exit(0);
  
}

unsigned int find_vec_location_double(std::vector<Double_t> vec, Double_t val)
{
  UInt_t size = vec.size();
  //for (int i=0;i<size;i++){cout << vec[i] << endl;}

  std::vector<Double_t>::iterator it = std::find(vec.begin(),vec.end(),val);
  if (it!=vec.end()) {
    return it - vec.begin();
  }
  else cout << "Can't locate " << val << " in your vector!" << endl; exit(0);
  
}

vector < Double_t > GetAlphaValues(Int_t runPeriod)
{
  Char_t temp[500];
  vector < Double_t > alphas (8,0.);
  sprintf(temp,"../smeared_peaks/nPE_Kev_%i.dat",runPeriod);
  ifstream infile;
  infile.open(temp);
  Int_t i = 0;
  
  while (infile >> alphas[i]) i++;
  return alphas;
}

void MB_calc_residuals_EvisPMTbyPMT(Int_t runPeriod)
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

  //Run Range for this envelope
  //Int_t runLow = 17359;
  //Int_t runHigh = 19959;
  Int_t calibrationPeriod = runPeriod;

  // Read data file
  char temp[500];
  sprintf(temp, "../residuals/source_runs_EvisPMTbyPMT_RunPeriod_%i.dat",calibrationPeriod);
  ifstream filein(temp);

  Int_t i = 0;
  Double_t run[500];
  string sourceName[500];
  Double_t EvisE1[500], EvisE2[500], EvisE3[500], EvisE4[500];
  Double_t EvisW1[500], EvisW2[500], EvisW3[500], EvisW4[500];
  Double_t resE1[500], resE2[500], resE3[500], resE4[500];
  Double_t resW1[500], resW2[500], resW3[500], resW4[500];
  Double_t err[500];
  while (!filein.eof()) {
    filein >> run[i] >> sourceName[i]
           >> EvisE1[i] >> EvisE2[i] >> EvisE3[i] >> EvisE4[i]
           >> EvisW1[i] >> EvisW2[i] >> EvisW3[i] >> EvisW4[i];
    if (filein.fail()) break;
    i++;
  }
  Int_t num = i;
  cout << "Number of data points: " << num << endl;

  // Load the smeared EQ values which are different for each PMT and source
  vector < vector <double> > EQsmeared = returnPeaks(calibrationPeriod,"EQ");
  for (int m=0; m<EQsmeared.size(); m++) {
    for (int mm=0; mm<EQsmeared[m].size(); mm++) {
      cout << EQsmeared[m][mm] << " ";
    }
    cout << endl;
  }
 
  //Load the simulated relationship between EQ and Etrue
  vector < vector <double> > EQ2Etrue = EQ2EtrueFit(calibrationPeriod);
  for (int m=0; m<EQ2Etrue.size(); m++) {
    for (int mm=0; mm<EQ2Etrue[m].size(); mm++) {
      cout << EQ2Etrue[m][mm] << " ";
    }
    cout << endl;
  }
  //Read in PMT quality file to save whether or not to use a PMT for each run
  vector<vector<int> > pmtQuality;
  vector <int> pmthold(8,0);
  vector<Int_t> pmtRun;
  sprintf(temp,"../residuals/PMT_runQuality_SrcPeriod_%i.dat",calibrationPeriod);
  ifstream pmt;
  pmt.open(temp);
  Int_t run_hold, EPMT1,EPMT2,EPMT3,EPMT4,WPMT1,WPMT2,WPMT3,WPMT4;
  Int_t numRuns=0;
  while (pmt >> run_hold >> pmthold[0] >> pmthold[1] >> pmthold[2]
	 >> pmthold[3] >> pmthold[4] >> pmthold[5]
	 >> pmthold[6] >> pmthold[7]) {
    pmtRun.push_back(run_hold);
    pmtQuality.push_back(pmthold);

//>> EPMT1 >> EPMT2 >> EPMT3 >> EPMT4 >> WPMT1 >> WPMT2 >> WPMT3 >> WPMT4;
    //cout << run_hold << endl;
    
    numRuns++;
    if (pmt.fail()) break;
  }
  pmt.close();
  for (int i=0;i<pmtRun.size();i++) {
    cout << pmtRun[i] << " " << pmtQuality[i][0] << " " << pmtQuality[i][1] << " " << pmtQuality[i][2] << 
      " " << pmtQuality[i][3] << " " << pmtQuality[i][4] << " " << pmtQuality[i][5] << " " << pmtQuality[i][6] << 
      " " << pmtQuality[i][7] << endl;
  }

  // Filling vectors with the weights for each PMT to be used when calculating the weighted average energy
  vector < Double_t > alphas = GetAlphaValues(calibrationPeriod);

  vector < vector<Double_t> > nPE_per_channel;
  Double_t nPE_hold, sigma_hold, mean_hold, nPE_per_chan_hold;
  nPE_per_channel.resize(pmtRun.size(), vector <Double_t> (8,0.));
  ifstream weightFile;
  for (int j=0; j<pmtRun.size(); j++) {
    sprintf(temp,"%s/nPE_weights_%i.dat",getenv("NPE_WEIGHTS"),pmtRun[j]);
    weightFile.open(temp);
    int i=0;
    while (weightFile >> mean_hold >> sigma_hold >> nPE_hold >> nPE_per_chan_hold) {
      nPE_per_channel[j][i]=nPE_per_chan_hold;
      i++;
      if (weightFile.fail()) break;
    }
    weightFile.close();
  }

  //cout << pmtRun[5] << endl;
  //for (int i=0;i<8;i++) {
  //  cout << nPE_per_channel[5][i] << endl;
  // }

  // Filling new vectors for each PMT with EVIS values only when the PMT is 
  // usable
  vector<int> runE1,runE2,runE3,runE4,runW1,runW2,runW3,runW4;
  vector<Double_t> EQE1,EQE2,EQE3,EQE4,EQW1,EQW2,EQW3,EQW4;
  vector<string> nameE1,nameE2,nameE3,nameE4,nameW1,nameW2,nameW3,nameW4;
  vector<Double_t> EVISE1, EVISE2, EVISE3, EVISE4;
  vector<Double_t> EVISW1, EVISW2, EVISW3, EVISW4;
  vector<Double_t> ResE1, ResE2, ResE3, ResE4;
  vector<Double_t> ResW1, ResW2, ResW3, ResW4;
  // Fill run, Evis, and eQ vectors
  for (Int_t i=0; i<num;i++) {
    UInt_t runPos = find_vec_location_int(pmtRun,(int)run[i]);
    cout << "Found run " << (int) run[i] << " in PMT Quality\n";
    int src_hold;
    if (sourceName[i]=="Ce") src_hold = 0;
    else if (sourceName[i]=="Sn") src_hold=1;
    else if (sourceName[i]=="Bi2") src_hold=2;
    else if (sourceName[i]=="Bi1") src_hold=3;

    if (pmtQuality[runPos][0]) {
      runE1.push_back((int)run[i]);
      EQE1.push_back(EQsmeared[0][src_hold]);
      nameE1.push_back(sourceName[i]);
      EVISE1.push_back(EvisE1[i]);
    }
    if (pmtQuality[runPos][1]) {
      runE2.push_back((int)run[i]);
      EQE2.push_back(EQsmeared[1][src_hold]);
      nameE2.push_back(sourceName[i]);
      EVISE2.push_back(EvisE2[i]);
    }
    if (pmtQuality[runPos][2]) {
      runE3.push_back((int)run[i]);
      EQE3.push_back(EQsmeared[2][src_hold]);
      nameE3.push_back(sourceName[i]);
      EVISE3.push_back(EvisE3[i]);
    }
    if (pmtQuality[runPos][3]) {
      runE4.push_back((int)run[i]);
      EQE4.push_back(EQsmeared[3][src_hold]);
      nameE4.push_back(sourceName[i]);
      EVISE4.push_back(EvisE4[i]);
    }
    if (pmtQuality[runPos][4]) {
      runW1.push_back((int)run[i]);
      EQW1.push_back(EQsmeared[4][src_hold]);
      nameW1.push_back(sourceName[i]);
      EVISW1.push_back(EvisW1[i]);
    }
    if (pmtQuality[runPos][5]) {
      runW2.push_back((int)run[i]);
      EQW2.push_back(EQsmeared[5][src_hold]);
      nameW2.push_back(sourceName[i]);
      EVISW2.push_back(EvisW2[i]);
    }
    if (pmtQuality[runPos][6]) {
      runW3.push_back((int)run[i]);
      EQW3.push_back(EQsmeared[6][src_hold]);
      nameW3.push_back(sourceName[i]);
      EVISW3.push_back(EvisW3[i]);
    }
    if (pmtQuality[runPos][7]) {
      runW4.push_back((int)run[i]);
      EQW4.push_back(EQsmeared[7][src_hold]);
      nameW4.push_back(sourceName[i]);
      EVISW4.push_back(EvisW4[i]);
    }
  }


  //Resize res vectors
  ResE1.resize(runE1.size());
  ResE2.resize(runE2.size());
  ResE3.resize(runE3.size());
  ResE4.resize(runE4.size());
  ResW1.resize(runW1.size());
  ResW2.resize(runW2.size());
  ResW3.resize(runW3.size());
  ResW4.resize(runW4.size());
  
  //Checking that there are enough sources to make a linearity curve
  /*if (std::find(EQE1.begin(),EQE1.end(),peakCe)==EQE1.end() || std::find(EQE1.begin(),EQE1.end(),peakSn)==EQE1.end() || std::find(EQE1.begin(),EQE1.end(),peakBiHigh)==EQE1.end()) { 
    cout << "Not enough sources to construct quadratic linearity curve\n";
    num=0;
    }*/
  if (std::find(nameE1.begin(),nameE1.end(),"Ce")==nameE1.end() || std::find(nameE1.begin(),nameE1.end(),"Sn")==nameE1.end() || std::find(nameE1.begin(),nameE1.end(),"Bi1")==nameE1.end()) { 
    cout << "Not enough sources to construct quadratic linearity curve\n";
    num=0;
  }
 

  // East 1

  sprintf(temp,"../residuals/residuals_EvisPMTbyPMT_East_runPeriod_%i_PMTE1.dat",calibrationPeriod);
  ofstream oFileE1(temp);

  if (runE1.size()>0 && (std::find(nameE1.begin(),nameE1.end(),"Ce")!=nameE1.end() && std::find(nameE1.begin(),nameE1.end(),"Sn")!=nameE1.end() && std::find(nameE1.begin(),nameE1.end(),"Bi1")!=nameE1.end())) {
 
    

    // Calculate residuals in [keV]
    
    for (int j=0; j<runE1.size(); j++) {
      if (nameE1[j]=="Ce") {
	ResE1[j] = EVISE1[j] - EQE1[j];
	oFileE1 << "Ce_East" << " " << runE1[j] << " " << ResE1[j] << endl;
	if (ResE1[j]>0.05*EQE1[j]) cout << "Ce_East" << " " << runE1[j] << " " << ResE1[j] << endl;
	//ResE1[j] = (fitEQ - peakCe)/peakCe * 100.;
      }
      else if (nameE1[j]=="Sn") {
	ResE1[j] = EVISE1[j] - EQE1[j];
	//ResE1[j] = (fitEQ - peakSn)/peakSn * 100.;
	oFileE1 << "Sn_East" << " " << runE1[j] << " " << ResE1[j] << endl;
	if (ResE1[j]>0.05*EQE1[j]) cout << "Sn_East" << " " << runE1[j] << " " << ResE1[j] << endl;
      }
      else if (nameE1[j]=="Bi1") {
	ResE1[j] = EVISE1[j] - EQE1[j];
	//ResE1[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE1 << "Bi1_East" << " " << runE1[j] << " " << ResE1[j] << endl;
	if (ResE1[j]>0.05*EQE1[j]) cout << "Bi1_East" << " " << runE1[j] << " " << ResE1[j] << endl;
 
      }
      else if (nameE1[j]=="Bi2") {
	ResE1[j] = EVISE1[j] - EQE1[j];
	//ResE1[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE1 << "Bi2_East" << " " << runE1[j] << " " << ResE1[j] << endl;
	if (ResE1[j]>0.05*EQE1[j]) cout << "Bi2_East" << " " << runE1[j] << " " << ResE1[j] << endl;
 
      }
    }

    // East 1 Residuals
    c1r = new TCanvas("c1r", "c1r");
    c1r->SetLogy(0);

    TGraphErrors *grE1r = new TGraphErrors(runE1.size(),&EVISE1[0],&ResE1[0],err,err);
    grE1r->SetTitle("");
    grE1r->SetMarkerColor(1);
    grE1r->SetLineColor(1);
    grE1r->SetMarkerStyle(20);
    grE1r->SetMarkerSize(0.75);
    grE1r->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
    grE1r->GetXaxis()->SetTitleOffset(1.2);
    grE1r->GetXaxis()->CenterTitle();
    grE1r->GetYaxis()->SetTitle("Residuals [keV]");
    grE1r->GetYaxis()->SetTitleOffset(1.2);
    grE1r->GetYaxis()->CenterTitle();
    grE1r->GetXaxis()->SetLimits(0.0,2400.0);
    grE1r->SetMinimum(-50.0);
    grE1r->SetMaximum( 50.0);
    grE1r->Draw("AP");

    const Int_t n = 2;
    Double_t x[n] = {0, 2400};
    Double_t y[n] = {0.0, 0.0};

    TGraph *gr1 = new TGraph(n,x,y);
    gr1->Draw("Same");
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->SetLineStyle(2);

    Double_t x1_text = 1200;
    Double_t y1_text = -40;

    TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("East PMT 1");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();
  }

  else {
    oFileE1 << "PMT NOT USABLE";}
  oFileE1.close();

  // East 2

  sprintf(temp,"../residuals/residuals_EvisPMTbyPMT_East_runPeriod_%i_PMTE2.dat",calibrationPeriod);
  ofstream oFileE2(temp);
  vector <Double_t> fitEQ_E2(runE2.size(),0);

  if (runE2.size()>0 && (std::find(nameE2.begin(),nameE2.end(),"Ce")!=nameE2.end() && std::find(nameE2.begin(),nameE2.end(),"Sn")!=nameE2.end() && std::find(nameE2.begin(),nameE2.end(),"Bi1")!=nameE2.end())) {

    

    // Calculate residuals in [keV]
    
    for (int j=0; j<runE2.size(); j++) {
 
      if (nameE2[j]=="Ce") {
	ResE2[j] = EVISE2[j] - EQE2[j];
	//ResE2[j] = (fitEQ - peakCe)/peakCe * 100.;
	oFileE2 << "Ce_East" << " " << runE2[j] << " " << ResE2[j] << endl;
	if (ResE2[j]>0.05*EQE2[j]) cout<< "Ce_East" << " " << runE2[j] << " " << ResE2[j] << endl;  
      }
      else if (nameE2[j]=="Sn") {
	ResE2[j] = EVISE2[j] - EQE2[j];
	//ResE2[j] = (fitEQ - peakSn)/peakSn * 100.;
	oFileE2 << "Sn_East" << " " << runE2[j] << " " << ResE2[j] << endl;
	if (ResE2[j]>0.05*EQE2[j]) cout<< "Sn_East" << " " << runE2[j] << " " << ResE2[j] << endl;  
      }
      else if (nameE2[j]=="Bi1") {
	ResE2[j] = EVISE2[j] - EQE2[j];
	//ResE2[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE2 << "Bi1_East" << " " << runE2[j] << " " << ResE2[j] << endl;
	if (ResE2[j]>0.05*EQE2[j]) cout<< "Bi1_East" << " " << runE2[j] << " " << ResE2[j] << endl;  
      }
      else if (nameE2[j]=="Bi2") {
	ResE2[j] = EVISE2[j] - EQE2[j];
	//ResE2[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE2 << "Bi2_East" << " " << runE2[j] << " " << ResE2[j] << endl;
	if (ResE2[j]>0.05*EQE2[j]) cout<< "Bi2_East" << " " << runE2[j] << " " << ResE2[j] << endl;  
      }
    }

    // East 2 residuals
    c2r = new TCanvas("c2r", "c2r");
    c2r->SetLogy(0);

    TGraphErrors *grE2r = new TGraphErrors(runE2.size(),&EVISE2[0],&ResE2[0],err,err);
    grE2r->SetTitle("");
    grE2r->SetMarkerColor(1);
    grE2r->SetLineColor(1);
    grE2r->SetMarkerStyle(20);
    grE2r->SetMarkerSize(0.75);
    grE2r->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
    grE2r->GetXaxis()->SetTitleOffset(1.2);
    grE2r->GetXaxis()->CenterTitle();
    grE2r->GetYaxis()->SetTitle("Residuals [keV]");
    grE2r->GetYaxis()->SetTitleOffset(1.2);
    grE2r->GetYaxis()->CenterTitle();
    grE2r->GetXaxis()->SetLimits(0.0,2400.0);
    grE2r->SetMinimum(-50.0);
    grE2r->SetMaximum( 50.0);
    grE2r->Draw("AP");

    Int_t n = 2;
    Double_t x[n] = {0, 2400};
    Double_t y[n] = {0.0, 0.0};

    TGraph *gr1 = new TGraph(n,x,y);
    gr1->Draw("Same");
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->SetLineStyle(2);

    Double_t x1_text = 1200;
    Double_t y1_text = -40;

    TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("East PMT 2");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();
  }

  else { 
    oFileE2 << "PMT NOT USABLE";}
  oFileE2.close();

  // East 3

  sprintf(temp,"../residuals/residuals_EvisPMTbyPMT_East_runPeriod_%i_PMTE3.dat",calibrationPeriod);
  ofstream oFileE3(temp);
  vector <Double_t> fitEQ_E3(runE3.size(),0);

  if (runE3.size()>0 && (std::find(nameE3.begin(),nameE3.end(),"Ce")!=nameE3.end() && std::find(nameE3.begin(),nameE3.end(),"Sn")!=nameE3.end() && std::find(nameE3.begin(),nameE3.end(),"Bi1")!=nameE3.end())) {

    

    // Calculate residuals in [keV]
    
    for (int j=0; j<runE3.size(); j++) {
    
      if (nameE3[j]=="Ce") {
	ResE3[j] = EVISE3[j] - EQE3[j];
	//ResE3[j] = (fitEQ - peakCe)/peakCe * 100.;
	oFileE3 << "Ce_East" << " " << runE3[j] << " " << ResE3[j] << endl;
	if (ResE3[j]>0.05*EQE3[j]) cout<< "Ce_East" << " " << runE3[j] << " " << ResE3[j] << endl;  
      }
      else if (nameE3[j]=="Sn") {
	ResE3[j] = EVISE3[j] - EQE3[j];
	//ResE3[j] = (fitEQ - peakSn)/peakSn * 100.;
	oFileE3 << "Sn_East" << " " << runE3[j] << " " << ResE3[j] << endl;
	if (ResE3[j]>0.05*EQE3[j]) cout << "Sn_East" << " " << runE3[j] << " " << ResE3[j] << endl;  
      }
      else if (nameE3[j]=="Bi1") {
	ResE3[j] = EVISE3[j] - EQE3[j];
	//ResE3[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE3 << "Bi1_East" << " " << runE3[j] << " " << ResE3[j] << endl;
	if (ResE3[j]>0.05*EQE3[j]) cout << "Bi1_East" << " " << runE3[j] << " " << ResE3[j] << endl;  
      }
      else if (nameE3[j]=="Bi2") {
	ResE3[j] = EVISE3[j] - EQE3[j];
	//ResE3[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE3 << "Bi2_East" << " " << runE3[j] << " " << ResE3[j] << endl;
	if (ResE3[j]>0.05*EQE3[j]) cout << "Bi2_East" << " " << runE3[j] << " " << ResE3[j] << endl;  
      }
    }

    // East 3 residuals
    c3r = new TCanvas("c3r", "c3r");
    c3r->SetLogy(0);

    TGraphErrors *grE3r = new TGraphErrors(runE3.size(),&EVISE3[0],&ResE3[0],err,err);
    grE3r->SetTitle("");
    grE3r->SetMarkerColor(1);
    grE3r->SetLineColor(1);
    grE3r->SetMarkerStyle(20);
    grE3r->SetMarkerSize(0.75);
    grE3r->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
    grE3r->GetXaxis()->SetTitleOffset(1.2);
    grE3r->GetXaxis()->CenterTitle();
    grE3r->GetYaxis()->SetTitle("Residuals [keV]");
    grE3r->GetYaxis()->SetTitleOffset(1.2);
    grE3r->GetYaxis()->CenterTitle();
    grE3r->GetXaxis()->SetLimits(0.0,2400.0);
    grE3r->SetMinimum(-50.0);
    grE3r->SetMaximum( 50.0);
    grE3r->Draw("AP");

    Int_t n = 2;
    Double_t x[n] = {0, 2400};
    Double_t y[n] = {0.0, 0.0};

    TGraph *gr1 = new TGraph(n,x,y);
    gr1->Draw("Same");
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->SetLineStyle(2);

    Double_t x1_text = 1200;
    Double_t y1_text = -40;

    TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("East PMT 3");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();
  }

  else {
    oFileE3 << "PMT NOT USABLE";}
  oFileE3.close();

  // East 4

  sprintf(temp,"../residuals/residuals_EvisPMTbyPMT_East_runPeriod_%i_PMTE4.dat",calibrationPeriod);
  ofstream oFileE4(temp);

  if (runE4.size()>0 && (std::find(nameE4.begin(),nameE4.end(),"Ce")!=nameE4.end() && std::find(nameE4.begin(),nameE4.end(),"Sn")!=nameE4.end() && std::find(nameE4.begin(),nameE4.end(),"Bi1")!=nameE4.end())) {

    
    // Calculate residuals in [keV]
    
    for (int j=0; j<runE4.size(); j++) {
      
      if (nameE4[j]=="Ce") {
	ResE4[j] = EVISE4[j] - EQE4[j];
	//ResE4[j] = (fitEQ - peakCe)/peakCe * 100.;
	oFileE4 << "Ce_East" << " " << runE4[j] << " " << ResE4[j] << endl;
	if (ResE4[j]>0.05*EQE4[j]) cout << "Ce_East" << " " << runE4[j] << " " << ResE4[j] << endl;  
      }
      else if (nameE4[j]=="Sn") {
	ResE4[j] = EVISE4[j] - EQE4[j];
	//ResE4[j] = (fitEQ - peakSn)/peakSn * 100.;
	oFileE4 << "Sn_East" << " " << runE4[j] << " " << ResE4[j] << endl;
	if (ResE4[j]>0.05*EQE4[j]) cout<< "Sn_East" << " " << runE4[j] << " " << ResE4[j] << endl;  
      }
      else if (nameE4[j]=="Bi1") {
	ResE4[j] = EVISE4[j] - EQE4[j];
	//ResE4[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE4 << "Bi1_East" << " " << runE4[j] << " " << ResE4[j] << endl;
	if (ResE4[j]>0.05*EQE4[j]) cout<< "Bi1_East" << " " << runE4[j] << " " << ResE4[j] << endl;  
      }
      else if (nameE4[j]=="Bi2") {
	ResE4[j] = EVISE4[j] - EQE4[j];
	//ResE4[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE4 << "Bi2_East" << " " << runE4[j] << " " << ResE4[j] << endl;
	if (ResE4[j]>0.05*EQE4[j]) cout<< "Bi2_East" << " " << runE4[j] << " " << ResE4[j] << endl;  
      }
    }

    // East 4 residuals
    c4r = new TCanvas("c4r", "c4r");
    c4r->SetLogy(0);

    TGraphErrors *grE4r = new TGraphErrors(runE4.size(),&EVISE4[0],&ResE4[0],err,err);
    grE4r->SetTitle("");
    grE4r->SetMarkerColor(1);
    grE4r->SetLineColor(1);
    grE4r->SetMarkerStyle(20);
    grE4r->SetMarkerSize(0.75);
    grE4r->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
    grE4r->GetXaxis()->SetTitleOffset(1.2);
    grE4r->GetXaxis()->CenterTitle();
    grE4r->GetYaxis()->SetTitle("Residuals [keV]");
    grE4r->GetYaxis()->SetTitleOffset(1.2);
    grE4r->GetYaxis()->CenterTitle();
    grE4r->GetXaxis()->SetLimits(0.0,2400.0);
    grE4r->SetMinimum(-50.0);
    grE4r->SetMaximum( 50.0);
    grE4r->Draw("AP");

    Int_t n = 2;
    Double_t x[n] = {0, 2400};
    Double_t y[n] = {0.0, 0.0};

    TGraph *gr1 = new TGraph(n,x,y);
    gr1->Draw("Same");
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->SetLineStyle(2);

    Double_t x1_text = 1200;
    Double_t y1_text = -40;

    TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("East PMT 4");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();
  }

  else {  oFileE4 << "PMT NOT USABLE";}
  oFileE4.close();

  ///////////////////////////////////////////////////////////////////////
  // Calcuting the weighted mean of the total energy of the East side  //
  ///////////////////////////////////////////////////////////////////////

  sprintf(temp,"../residuals/residuals_EvisPMTbyPMT_East_runPeriod_%i.dat",calibrationPeriod);
  ofstream oFileE(temp);

  cout << "CALCULATING EAST RESIDUALS" << endl;
  vector <Double_t> Etrue_East(num,0);
  vector <Double_t> x_East(num,0);
  vector <Double_t> res_East(num,0);

  std::vector<Int_t>::iterator E1 = runE1.begin();
  std::vector<Double_t>::iterator E1_EQ = EQE1.begin();
  std::vector<Double_t>::iterator E1_EVIS = EVISE1.begin();
  std::vector<Int_t>::iterator E2 = runE2.begin();
  std::vector<Double_t>::iterator E2_EQ = EQE2.begin();
  std::vector<Double_t>::iterator E2_EVIS = EVISE2.begin();
  std::vector<Int_t>::iterator E3 = runE3.begin();
  std::vector<Double_t>::iterator E3_EQ = EQE3.begin();
  std::vector<Double_t>::iterator E3_EVIS = EVISE3.begin();
  std::vector<Int_t>::iterator E4 = runE4.begin();
  std::vector<Double_t>::iterator E4_EQ = EQE4.begin();
  std::vector<Double_t>::iterator E4_EVIS = EVISE4.begin();
  std::vector<string>::iterator nmE1 = nameE1.begin();
  std::vector<string>::iterator nmE2 = nameE2.begin();
  std::vector<string>::iterator nmE3 = nameE3.begin();
  std::vector<string>::iterator nmE4 = nameE4.begin();
  std::vector<Double_t>::iterator RE1 = ResE1.begin();
  std::vector<Double_t>::iterator RE2 = ResE2.begin();
  std::vector<Double_t>::iterator RE3 = ResE3.begin();
  std::vector<Double_t>::iterator RE4 = ResE4.begin();

  
  for (int j=0;j<num;j++) {
    Double_t weight1, weight2, weight3, weight4; //these will hold the 4 weights
    Double_t Energy1, Energy2, Energy3, Energy4; //These will hold the energies for each pmt in true energy
    Double_t ResidualE1, ResidualE2, ResidualE3, ResidualE4; //These hold the residuals for the weighting method..
    UInt_t pmtVecLocation = find_vec_location_int(pmtRun,(int)run[j]);

    if (pmtQuality[pmtVecLocation][0] && *E1==(int)run[j] && *nmE1==sourceName[j]) {
      Double_t N = EvisE1[j]*alphas[0];
      Double_t f = sqrt(N)/N;
      ResidualE1 = *RE1;
      //Energy1 = EQ2Etrue[0][0]+EQ2Etrue[0][1]*(*E1_fitEQ)+EQ2Etrue[0][2]*(*E1_fitEQ)*(*E1_fitEQ);
      Energy1 = E1_EVIS;
      weight1 = 1/(Energy1*Energy1*f*f);
      std::advance(E1,1);
      std::advance(E1_EQ,1);
      std::advance(E1_EVIS,1);
      std::advance(nmE1,1);
      std::advance(RE1,1);
    }
    else {weight1=0.; Energy1=0.; ResidualE1=0.;}

    if (pmtQuality[pmtVecLocation][1] && *E2==(int)run[j] && *nmE2==sourceName[j]) {
      Double_t N = EvisE2[j]*alphas[1];
      Double_t f = sqrt(N)/N;
      ResidualE2 = *RE2;
      //Energy2 = EQ2Etrue[1][0]+EQ2Etrue[1][1]*(*E2_fitEQ)+EQ2Etrue[1][2]*(*E2_fitEQ)*(*E2_fitEQ);
      Energy2 = E2_EVIS;
      weight2 = 1/(Energy2*Energy2*f*f);
      std::advance(E2,1);
      std::advance(E2_EQ,1);
      std::advance(E2_EVIS,1);
      std::advance(nmE2,1);
      std::advance(RE2,1);
    }
    else {weight2=0.; Energy2=0.; ResidualE2=0.;}

    if (pmtQuality[pmtVecLocation][2] && *E3==(int)run[j] && *nmE3==sourceName[j]) {
      Double_t N = EvisE3[j]*alphas[2];
      Double_t f = sqrt(N)/N;
      ResidualE4 = *RE3;
      //Energy3 = EQ2Etrue[2][0]+EQ2Etrue[2][1]*(*E3_fitEQ)+EQ2Etrue[2][2]*(*E3_fitEQ)*(*E3_fitEQ);;
      Energy3 = E3_EVIS;
      weight3 = 1/(Energy3*Energy3*f*f);
      std::advance(E3,1);
      std::advance(E3_EQ,1);
      std::advance(E3_EVIS,1);
      std::advance(nmE3,1);
      std::advance(RE3,1);
    }
    else {weight3=0.; Energy3=0.; ResidualE3=0.;}

    if (pmtQuality[pmtVecLocation][3] && *E4==(int)run[j] && *nmE4==sourceName[j]) {
      Double_t N = EvisE4[j]*alphas[3];
      Double_t f = sqrt(N)/N;
      ResidualE4 = *RE4;
      //Energy4 = EQ2Etrue[3][0]+EQ2Etrue[3][1]*(*E4_fitEQ)+EQ2Etrue[3][2]*(*E4_fitEQ)*(*E4_fitEQ);
      Energy4 = E4_EVIS;
      weight4 = 1/(Energy4*Energy4*f*f);
      std::advance(E4,1);
      std::advance(E4_EQ,1);
      std::advance(E4_EVIS,1);
      std::advance(nmE4,1);
      std::advance(RE4,1);
    }
    else {weight4=0.; Energy4=0.; ResidualE4=0.;}
  
    Etrue_East[j] = (weight1*Energy1+weight2*Energy2+weight3*Energy3+weight4*Energy4)/(weight1+weight2+weight3+weight4);
    res_East[j] = (weight1*ResidualE1+weight2*ResidualE2+weight3*ResidualE3+weight4*ResidualE4)/(weight1+weight2+weight3+weight4);
    if (sourceName[j]=="Ce") {
      //res_East[j] = Etrue_East[j] - peakCe;
      x_East[j] = peakCe_EQ;
      oFileE << "Ce_East" << " " << (int) run[j] << " " << res_East[j] << endl;
    }
    else if (sourceName[j]=="Sn") {
      //res_East[j] = Etrue_East[j] - peakSn;
      x_East[j] = peakSn_EQ;
      oFileE << "Sn_East" << " " << (int) run[j] << " " << res_East[j] << endl;
      /*cout << Energy1 << " " << weight1 << " "
	   << Energy2 << " " << weight2 << " "
	   << Energy3 << " " << weight3 << " " 
	   << Energy4 << " " << weight4 << endl;*/
    }
    else if (sourceName[j]=="Bi1") {
      //res_East[j] = Etrue_East[j] - peakBiHigh;
      x_East[j] = peakBiHigh_EQ;
      oFileE << "Bi1_East" << " " << (int) run[j] << " " << res_East[j] << endl;
   }
    else if (sourceName[j]=="Bi2") {
      //res_East[j] = Etrue_East[j] - peakBiHigh;
      x_East[j] = peakBiLow_EQ;
      oFileE << "Bi2_East" << " " << (int) run[j] << " " << res_East[j] << endl;
   }

  }
  oFileE.close();
  
  // East Average residuals
  cEr = new TCanvas("cEr", "cEr");
  cEr->SetLogy(0);
  TGraphErrors *grEr; 

  if (num>0) {
    grEr = new TGraphErrors(num,&x_East[0],&res_East[0],err,err);
    grEr->SetTitle("");
    grEr->SetMarkerColor(1);
    grEr->SetLineColor(1);
    grEr->SetMarkerStyle(20);
    grEr->SetMarkerSize(0.75);
    grEr->GetXaxis()->SetTitle("E_{Q} [keV]");
    grEr->GetXaxis()->SetTitleOffset(1.2);
    grEr->GetXaxis()->CenterTitle();
    grEr->GetYaxis()->SetTitle("Residuals [keV]");
    grEr->GetYaxis()->SetTitleOffset(1.2);
    grEr->GetYaxis()->CenterTitle();
    grEr->GetXaxis()->SetLimits(0.0,1200.0);
    grEr->SetMinimum(-30.0);
    grEr->SetMaximum( 30.0);
    grEr->Draw("AP");
    
    Int_t n = 2;
    Double_t x[n] = {0, 2400};
    Double_t y[n] = {0.0, 0.0};
    
    TGraph *gr1 = new TGraph(n,x,y);
    gr1->Draw("Same");
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->SetLineStyle(2);
    
    Double_t x1_text = 1000;
    Double_t y1_text = -20;
    
    TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("East");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();
  }
  ////////////////////////////////////////////////////////////////////

  // West 1

  sprintf(temp,"../residuals/residuals_EvisPMTbyPMT_West_runPeriod_%i_PMTW1.dat",calibrationPeriod);
  ofstream oFileW1(temp);

  if (runW1.size()>0 && (std::find(nameW1.begin(),nameW1.end(),"Ce")!=nameW1.end() && std::find(nameW1.begin(),nameW1.end(),"Sn")!=nameW1.end() && std::find(nameW1.begin(),nameW1.end(),"Bi1")!=nameW1.end())) {

    
  
    // Calculate residuals in [keV]
    //Double_t fitEQ_W1[num];
    for (int j=0; j<runW1.size(); j++) {
      
      if (nameW1[j]=="Ce") {
	ResW1[j] = EVISW1[j] - EQW1[j];
	oFileW1 << "Ce_West" << " " << runW1[j] << " " << ResW1[j] << endl;
	if (ResW1[j]>0.05*EQW1[j]) cout<< "Ce_West" << " " << runW1[j] << " " << ResW1[j] << endl;  
      }
      else if (nameW1[j]=="Sn") {
	ResW1[j] = EVISW1[j] - EQW1[j];
	oFileW1 << "Sn_West" << " " <<  runW1[j] << " " << ResW1[j] << endl;
	if (ResW1[j]>0.05*EQW1[j]) cout<< "Sn_West" << " " << runW1[j] << " " << ResW1[j] << endl;  
      }
      else if (nameW1[j]=="Bi1") {
	ResW1[j] = EVISW1[j] - EQW1[j];
	oFileW1 << "Bi1_West" << " " << runW1[j] << " " << ResW1[j] << endl;
	if (ResW1[j]>0.05*EQW1[j]) cout<< "Bi1_West" << " " << runW1[j] << " " << ResW1[j] << endl;  
      }
      else if (nameW1[j]=="Bi2") {
	ResW1[j] = EVISW1[j] - EQW1[j];
	oFileW1 << "Bi2_West" << " " << runW1[j] << " " << ResW1[j] << endl;
	if (ResW1[j]>0.05*EQW1[j]) cout<< "Bi2_West" << " " << runW1[j] << " " << ResW1[j] << endl;  
      }
    }

    // West 1 residuals
    cW1r = new TCanvas("cW1r", "cW1r");
    cW1r->SetLogy(0);

    TGraphErrors *grW1r = new TGraphErrors(runW1.size(),&EVISW1[0],&ResW1[0],err,err);
    grW1r->SetTitle("");
    grW1r->SetMarkerColor(1);
    grW1r->SetLineColor(1);
    grW1r->SetMarkerStyle(20);
    grW1r->SetMarkerSize(0.75);
    grW1r->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
    grW1r->GetXaxis()->SetTitleOffset(1.2);
    grW1r->GetXaxis()->CenterTitle();
    grW1r->GetYaxis()->SetTitle("Residuals [keV]");
    grW1r->GetYaxis()->SetTitleOffset(1.2);
    grW1r->GetYaxis()->CenterTitle();
    grW1r->GetXaxis()->SetLimits(0.0,2400.0);
    grW1r->SetMinimum(-50.0);
    grW1r->SetMaximum( 50.0);
    grW1r->Draw("AP");

    Int_t n = 2;
    Double_t x[n] = {0, 2400};
    Double_t y[n] = {0.0, 0.0};

    TGraph *gr1 = new TGraph(n,x,y);
    gr1->Draw("Same");
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->SetLineStyle(2);

    Double_t x1_text = 1200;
    Double_t y1_text = -40;

    TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("West PMT 1");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();
  }

  else {
    oFileW1 << "PMT NOT USABLE";}
  oFileW1.close();

  // West 2

  sprintf(temp,"../residuals/residuals_EvisPMTbyPMT_West_runPeriod_%i_PMTW2.dat",calibrationPeriod);
  ofstream oFileW2(temp);
 

  if (runW2.size()>0 && (std::find(nameW2.begin(),nameW2.end(),"Ce")!=nameW2.end() && std::find(nameW2.begin(),nameW2.end(),"Sn")!=nameW2.end() && std::find(nameW2.begin(),nameW2.end(),"Bi1")!=nameW2.end())) {

  

    // Calculate residuals in [keV]
    //Double_t fitEQ_W2[num];
    for (int j=0; j<runW2.size(); j++) {
     
      if (nameW2[j]=="Ce") {
	ResW2[j] = EVISW2[j] - EQW2[j];
	oFileW2 << "Ce_West" << " " << runW2[j] << " " << ResW2[j] << endl;
	if (ResW2[j]>0.05*EQW2[j]) cout<< "Ce_West" << " " << runW2[j] << " " << ResW2[j] << endl;  
      }
      else if (nameW2[j]=="Sn") {
	ResW2[j] = EVISW2[j] - EQW2[j];
	oFileW2 << "Sn_West" << " " << runW2[j] << " " << ResW2[j] << endl;
	if (ResW2[j]>0.05*EQW2[j]) cout<< "Sn_West" << " " << runW2[j] << " " << ResW2[j] << endl;  
      }
      else if (nameW2[j]=="Bi1") {
	ResW2[j] = EVISW2[j] - EQW2[j];
	oFileW2 << "Bi1_West" << " " << runW2[j] << " " << ResW2[j] << endl;
	if (ResW2[j]>0.05*EQW2[j]) cout<< "Bi1_West" << " " << runW2[j] << " " << ResW2[j] << endl;  
      }
      else if (nameW2[j]=="Bi2") {
	ResW2[j] = EVISW2[j] - EQW2[j];
	oFileW2 << "Bi2_West" << " " << runW2[j] << " " << ResW2[j] << endl;
	if (ResW2[j]>0.05*EQW2[j]) cout<< "Bi2_West" << " " << runW2[j] << " " << ResW2[j] << endl;  
      }
    }

    // West 2 residuals
    cW2r = new TCanvas("cW2r", "cW2r");
    cW2r->SetLogy(0);

    TGraphErrors *grW2r = new TGraphErrors(runW2.size(),&EVISW2[0],&ResW2[0],err,err);
    grW2r->SetTitle("");
    grW2r->SetMarkerColor(1);
    grW2r->SetLineColor(1);
    grW2r->SetMarkerStyle(20);
    grW2r->SetMarkerSize(0.75);
    grW2r->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointEVIS [channels]");
    grW2r->GetXaxis()->SetTitleOffset(1.2);
    grW2r->GetXaxis()->CenterTitle();
    grW2r->GetYaxis()->SetTitle("Residuals [keV]");
    grW2r->GetYaxis()->SetTitleOffset(1.2);
    grW2r->GetYaxis()->CenterTitle();
    grW2r->GetXaxis()->SetLimits(0.0,2400.0);
    grW2r->SetMinimum(-50.0);
    grW2r->SetMaximum( 50.0);
    grW2r->Draw("AP");

    Int_t n = 2;
    Double_t x[n] = {0, 2400};
    Double_t y[n] = {0.0, 0.0};

    TGraph *gr1 = new TGraph(n,x,y);
    gr1->Draw("Same");
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->SetLineStyle(2);

    Double_t x1_text = 1200;
    Double_t y1_text = -40;

    TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("West PMT 2");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();
  }

  else {
    oFileW2 << "PMT NOT USABLE";}
  oFileW2.close();

  // West 3

  sprintf(temp,"../residuals/residuals_EvisPMTbyPMT_West_runPeriod_%i_PMTW3.dat",calibrationPeriod);
  ofstream oFileW3(temp);
  

  if (runW3.size()>0 && (std::find(nameW3.begin(),nameW3.end(),"Ce")!=nameW3.end() && std::find(nameW3.begin(),nameW3.end(),"Sn")!=nameW3.end() && std::find(nameW3.begin(),nameW3.end(),"Bi1")!=nameW3.end())) {

    

    // Calculate residuals in [keV]
    //Double_t fitEQ_W3[num];
    for (int j=0; j<runW3.size(); j++) {
     
      if (nameW3[j]=="Ce") {
	ResW3[j] = EVISW3[j] - EQW3[j];
	oFileW3 << "Ce_West" << " " << runW3[j] << " " << ResW3[j] << endl;
	if (ResW3[j]>0.05*EQW3[j]) cout<< "Ce_West" << " " << runW3[j] << " " << ResW3[j] << endl;  
      }
      else if (nameW3[j]=="Sn") {
	ResW3[j] = EVISW3[j] - EQW3[j];
	oFileW3 << "Sn_West" << " " << runW3[j] << " " << ResW3[j] << endl;
	if (ResW3[j]>0.05*EQW3[j]) cout<< "Sn_West" << " " << runW3[j] << " " << ResW3[j] << endl;  
      }
      else if (nameW3[j]=="Bi1") {
	ResW3[j] = EVISW3[j] - EQW3[j];
	oFileW3 << "Bi1_West" << " " << runW3[j] << " " << ResW3[j] << endl;
	if (ResW3[j]>0.05*EQW3[j]) cout<< "Bi1_West" << " " << runW3[j] << " " << ResW3[j] << endl;  
      }
      else if (nameW3[j]=="Bi2") {
	ResW3[j] = EVISW3[j] - EQW3[j];
	oFileW3 << "Bi2_West" << " " << runW3[j] << " " << ResW3[j] << endl;
	if (ResW3[j]>0.05*EQW3[j]) cout<< "Bi2_West" << " " << runW3[j] << " " << ResW3[j] << endl;  
      }
    }

    // West 3 residuals
    cW3r = new TCanvas("cW3r", "cW3r");
    cW3r->SetLogy(0);

    TGraphErrors *grW3r = new TGraphErrors(runW3.size(),&EVISW3[0],&ResW3[0],err,err);
    grW3r->SetTitle("");
    grW3r->SetMarkerColor(1);
    grW3r->SetLineColor(1);
    grW3r->SetMarkerStyle(20);
    grW3r->SetMarkerSize(0.75);
    grW3r->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
    grW3r->GetXaxis()->SetTitleOffset(1.2);
    grW3r->GetXaxis()->CenterTitle();
    grW3r->GetYaxis()->SetTitle("Residuals [keV]");
    grW3r->GetYaxis()->SetTitleOffset(1.2);
    grW3r->GetYaxis()->CenterTitle();
    grW3r->GetXaxis()->SetLimits(0.0,2400.0);
    grW3r->SetMinimum(-50.0);
    grW3r->SetMaximum( 50.0);
    grW3r->Draw("AP");

    Int_t n = 2;
    Double_t x[n] = {0, 2400};
    Double_t y[n] = {0.0, 0.0};

    TGraph *gr1 = new TGraph(n,x,y);
    gr1->Draw("Same");
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->SetLineStyle(2);

    Double_t x1_text = 1200;
    Double_t y1_text = -40;

    TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("West PMT 3");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();
  }

  else { 
    oFileW3 << "PMT NOT USABLE";}
  oFileW3.close();

  // West 4

  sprintf(temp,"../residuals/residuals_EvisPMTbyPMT_West_runPeriod_%i_PMTW4.dat",calibrationPeriod);
  ofstream oFileW4(temp);

  if (runW4.size()>0 && (std::find(nameW4.begin(),nameW4.end(),"Ce")!=nameW4.end() && std::find(nameW4.begin(),nameW4.end(),"Sn")!=nameW4.end() && std::find(nameW4.begin(),nameW4.end(),"Bi1")!=nameW4.end())) {

    
    // Calculate residuals in [keV]
    //Double_t fitEQ_W4[num];
    for (int j=0; j<runW4.size(); j++) {
  
      if (nameW4[j]=="Ce") {
	ResW4[j] = EVISW4[j] - EQW4[j];
	oFileW4 << "Ce_West" << " " << runW4[j] << " " << ResW4[j] << endl;
	if (ResW4[j]>0.05*EQW4[j]) cout<< "Ce_West" << " " << runW4[j] << " " << ResW4[j] << endl;  
      }
      else if (nameW4[j]=="Sn") {
	ResW4[j] = EVISW4[j] - EQW4[j];
	oFileW4 << "Sn_West" << " " << runW4[j] << " " << ResW4[j] << endl;
	if (ResW4[j]>0.05*EQW4[j]) cout<< "Sn_West" << " " << runW4[j] << " " << ResW4[j] << endl;  
      }
      else if (nameW4[j]=="Bi1"){
	ResW4[j] = EVISW4[j] - EQW4[j];
	oFileW4 << "Bi1_West" << " " << runW4[j] << " " << ResW4[j] << endl;
	if (ResW4[j]>0.05*EQW4[j]) cout<< "Bi1_West" << " " << runW4[j] << " " << ResW4[j] << endl;  
      }
      else if (nameW4[j]=="Bi2"){
	ResW4[j] = EVISW4[j] - EQW4[j];
	oFileW4 << "Bi2_West" << " " << runW4[j] << " " << ResW4[j] << endl;
	if (ResW4[j]>0.05*EQW4[j]) cout<< "Bi2_West" << " " << runW4[j] << " " << ResW4[j] << endl;  
      }
    }

    // West 4 residuals
    cW4r = new TCanvas("cW4r", "cW4r");
    cW4r->SetLogy(0);

    TGraphErrors *grW4r = new TGraphErrors(runW4.size(),&EVISW4[0],&ResW4[0],err,err);
    grW4r->SetTitle("");
    grW4r->SetMarkerColor(1);
    grW4r->SetLineColor(1);
    grW4r->SetMarkerStyle(20);
    grW4r->SetMarkerSize(0.75);
    grW4r->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
    grW4r->GetXaxis()->SetTitleOffset(1.2);
    grW4r->GetXaxis()->CenterTitle();
    grW4r->GetYaxis()->SetTitle("Residuals [keV]");
    grW4r->GetYaxis()->SetTitleOffset(1.2);
    grW4r->GetYaxis()->CenterTitle();
    grW4r->GetXaxis()->SetLimits(0.0,4000.);
    grW4r->SetMinimum(-50.0);
    grW4r->SetMaximum( 50.0);
    grW4r->Draw("AP");

    Int_t n = 2;
    Double_t x[n] = {0, 4000.};
    Double_t y[n] = {0.0, 0.0};

    TGraph *gr1 = new TGraph(n,x,y);
    gr1->Draw("Same");
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->SetLineStyle(2);

    Double_t x1_text = 1200;
    Double_t y1_text = -40;

    TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("West PMT 4");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();
  }
  
  else { 
    oFileW4 << "PMT NOT USABLE";
  }
  oFileW4.close();


  ///////////////////////////////////////////////////////////////////////
  // Calculating the weighted mean of the total energy of the West side  //
  ///////////////////////////////////////////////////////////////////////

  sprintf(temp,"../residuals/residuals_EvisPMTbyPMT_West_runPeriod_%i.dat",calibrationPeriod);
  ofstream oFileW(temp);

  cout << "CALCULATING WEST RESIDUALS" << endl;
  vector <Double_t> Etrue_West(num,0);
  vector <Double_t> x_West(num,0);
  vector <Double_t> res_West(num,0);

  std::vector<Int_t>::iterator W1 = runW1.begin();
  std::vector<Double_t>::iterator W1_EQ = EQW1.begin();
  std::vector<Double_t>::iterator W1_EVIS = EVISW1.begin();
  std::vector<Int_t>::iterator W2 = runW2.begin();
  std::vector<Double_t>::iterator W2_EQ = EQW2.begin();
  std::vector<Double_t>::iterator W2_EVIS = EVISW2.begin();
  std::vector<Int_t>::iterator W3 = runW3.begin();
  std::vector<Double_t>::iterator W3_EQ = EQW3.begin();
  std::vector<Double_t>::iterator W3_EVIS = EVISW3.begin();
  std::vector<Int_t>::iterator W4 = runW4.begin();
  std::vector<Double_t>::iterator W4_EQ = EQW4.begin();
  std::vector<Double_t>::iterator W4_EVIS = EVISW4.begin();
  std::vector<string>::iterator nmW1 = nameW1.begin();
  std::vector<string>::iterator nmW2 = nameW2.begin();
  std::vector<string>::iterator nmW3 = nameW3.begin();
  std::vector<string>::iterator nmW4 = nameW4.begin();
  std::vector<Double_t>::iterator RW1 = ResW1.begin();
  std::vector<Double_t>::iterator RW2 = ResW2.begin();
  std::vector<Double_t>::iterator RW3 = ResW3.begin();
  std::vector<Double_t>::iterator RW4 = ResW4.begin();

  for (int j=0;j<num;j++) {
    Double_t weight1, weight2, weight3, weight4; //these will hold the 4 weights
    Double_t Energy1, Energy2, Energy3, Energy4; //These will hold the energies for each pmt
    Double_t ResidualW1, ResidualW2, ResidualW3, ResidualW4;
    UInt_t pmtVecLocation = find_vec_location_int(pmtRun,(int)run[j]);

    if (pmtQuality[pmtVecLocation][4] && *W1==(int)run[j] && *nmW1==sourceName[j]) {
      Double_t N = EvisW1[j]*alphas[4];
      Double_t f = sqrt(N)/N;
      ResidualW1 = *RW1;
      Energy1 = W1_EVIS;
      std::advance(RW1,1);
      //Energy1 = EQ2Etrue[4][0]+EQ2Etrue[4][1]*(*W1_fitEQ)+EQ2Etrue[4][2]*(*W1_EVIS)*(*W1_fitEQ);
      weight1 = 1/(Energy1*Energy1*f*f);
      std::advance(W1,1);
      std::advance(W1_EQ,1);
      std::advance(W1_EVIS,1);
      std::advance(nmW1,1);
    }
    else {weight1=0.; Energy1=0.; ResidualW1=0.;}

    if (pmtQuality[pmtVecLocation][5] && *W2==(int)run[j] && *nmW2==sourceName[j]) {
      Double_t N = EvisW2[j]*alphas[5];
      Double_t f = sqrt(N)/N;
      ResidualW2 = *RW2;
      Energy2 = W2_EVIS;
      std::advance(RW2,1);
      //Energy2 = EQ2Etrue[5][0]+EQ2Etrue[5][1]*(*W2_fitEQ)+EQ2Etrue[5][2]*(*W2_fitEQ)*(*W2_fitEQ);
      weight2 = 1/(Energy2*Energy2*f*f);
      std::advance(W2,1);
      std::advance(W2_EQ,1);
      std::advance(W2_EVIS,1);
      std::advance(nmW2,1);
    }
    else {weight2=0.; Energy2=0.; ResidualW2=0.;}

    if (pmtQuality[pmtVecLocation][6] && *W3==(int)run[j] && *nmW3==sourceName[j]) {
      Double_t N = EvisW3[j]*alphas[6];
      Double_t f = sqrt(N)/N;
      ResidualW3 = *RW3;
      Energy3 = W3_EVIS;
      std::advance(RW3,1);
      //Energy3 = EQ2Etrue[6][0]+EQ2Etrue[6][1]*(*W3_fitEQ)+EQ2Etrue[6][2]*(*W3_fitEQ)*(*W3_fitEQ);
      weight3 = 1/(Energy3*Energy3*f*f);
      std::advance(W3,1);
      std::advance(W3_EQ,1);
      std::advance(W3_EVIS,1);
      std::advance(nmW3,1);
    }
    else {weight3=0.; Energy3=0.; ResidualW3=0.;}

    if (pmtQuality[pmtVecLocation][7] && *W4==(int)run[j] && *nmW4==sourceName[j]) {
      Double_t N = EvisW4[j]*alphas[7];
      Double_t f = sqrt(N)/N;
      ResidualW4 = *RW4;
      Energy4 = W4_EVIS;
      std::advance(RW4,1);
      //Energy4 = EQ2Etrue[7][0]+EQ2Etrue[7][1]*(*W4_fitEQ)+EQ2Etrue[7][2]*(*W4_fitEQ)*(*W4_fitEQ);
      weight4 = 1/(Energy4*Energy4*f*f);
      std::advance(W4,1);
      std::advance(W4_EQ,1);
      std::advance(W4_EVIS,1);
      std::advance(nmW4,1);
    }
    else {weight4=0.; Energy4=0.; ResidualW4=0.;}

    Etrue_West[j] = (weight1*Energy1+weight2*Energy2+weight3*Energy3+weight4*Energy4)/(weight1+weight2+weight3+weight4);
    res_West[j] = (weight1*ResidualW1+weight2*ResidualW2+weight3*ResidualW3+weight4*ResidualW4)/(weight1+weight2+weight3+weight4);
    if (sourceName[j]=="Ce") {
      //res_West[j] = Etrue_West[j] - peakCe;
      x_West[j] = peakCe_EQ;
      oFileW << "Ce_West" << " " << (int) run[j] << " " << res_West[j] << endl;
    }
    else if (sourceName[j]=="Sn") {
      //res_West[j] = Etrue_West[j] - peakSn;
      x_West[j] = peakSn_EQ;
      oFileW << "Sn_West" << " " << (int) run[j] << " " << res_West[j] << endl;
      /*cout << Energy1 << " " << weight1 << " "
	   << Energy2 << " " << weight2 << " "
	   << Energy3 << " " << weight3 << " " 
	   << Energy4 << " " << weight4 << endl;*/
    }
    else if (sourceName[j]=="Bi1") {
      //res_West[j] = Etrue_West[j] - peakBiHigh;
      x_West[j] = peakBiHigh_EQ;
      oFileW << "Bi1_West" << " " << (int) run[j] << " " << res_West[j] << endl;
   }
    else if (sourceName[j]=="Bi2") {
      //res_West[j] = Etrue_West[j] - peakBiHigh;
      x_West[j] = peakBiLow_EQ;
      oFileW << "Bi2_West" << " " << (int) run[j] << " " << res_West[j] << endl;
   }

  }
  oFileW.close();

  // West Average residuals
  cWr = new TCanvas("cWr", "cWr");
  cWr->SetLogy(0);
  TGraphErrors *grWr;
 
  if (num>0) {
    grWr = new TGraphErrors(num,&x_West[0],&res_West[0],err,err);
    grWr->SetTitle("");
    grWr->SetMarkerColor(1);
    grWr->SetLineColor(1);
    grWr->SetMarkerStyle(20);
    grWr->SetMarkerSize(0.75);
    grWr->GetXaxis()->SetTitle("E_{Q} [keV]");
    grWr->GetXaxis()->SetTitleOffset(1.2);
    grWr->GetXaxis()->CenterTitle();
    grWr->GetYaxis()->SetTitle("Residuals [keV]");
    grWr->GetYaxis()->SetTitleOffset(1.2);
    grWr->GetYaxis()->CenterTitle();
    grWr->GetXaxis()->SetLimits(0.0,1200.0);
    grWr->SetMinimum(-30.0);
    grWr->SetMaximum( 30.0);
    grWr->Draw("AP");
    
    Int_t n = 2;
    Double_t x[n] = {0, 2400};
    Double_t y[n] = {0.0, 0.0};
    
    TGraph *gr1 = new TGraph(n,x,y);
    gr1->Draw("Same");
    gr1->SetLineWidth(2);
    gr1->SetLineColor(1);
    gr1->SetLineStyle(2);
    
    Double_t x1_text = 1000;
    Double_t y1_text = -20;
    
    TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("West");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();
  }
  
  //outfile->Write();
  //outfile->Close();

}
