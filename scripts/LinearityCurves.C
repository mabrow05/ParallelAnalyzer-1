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

void LinearityCurves(Int_t runPeriod)
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
  sprintf(temp, "../residuals/source_runs_RunPeriod_%i.dat",calibrationPeriod);
  ifstream filein(temp);
  sprintf(temp, "../residuals/SIM_source_runs_EvisPMTbyPMT_RunPeriod_%i.dat",calibrationPeriod);
  ifstream simFileIn(temp);

  Int_t i = 0;
  Int_t run[500], run_hold=0;
  string sourceName[500], sourceName_hold="";
  Double_t adcE1[500], adcE2[500], adcE3[500], adcE4[500];
  Double_t adcW1[500], adcW2[500], adcW3[500], adcW4[500];
  Double_t EqE1[500], EqE2[500], EqE3[500], EqE4[500];
  Double_t EqW1[500], EqW2[500], EqW3[500], EqW4[500];
  Double_t resE1[500], resE2[500], resE3[500], resE4[500];
  Double_t resW1[500], resW2[500], resW3[500], resW4[500];
  Double_t err[500];
  Double_t Eq_hold[8];
  
  while (filein >> run[i] >> sourceName[i]
	 >> adcE1[i] >> adcE2[i] >> adcE3[i] >> adcE4[i]
	 >> adcW1[i] >> adcW2[i] >> adcW3[i] >> adcW4[i]) {

    while (!(run_hold==run[i] && sourceName_hold==sourceName[i])) {
      simFileIn >> run_hold >> sourceName_hold >> Eq_hold[0] >> Eq_hold[1] >> Eq_hold[2]
		>> Eq_hold[3] >> Eq_hold[4] >> Eq_hold[5] >> Eq_hold[6] >> Eq_hold[7];
    }
    //cout << run_hold << " " << sourceName_hold << " ";
    if (run_hold==run[i] && sourceName_hold==sourceName[i]) { 
      EqE1[i] = Eq_hold[0];
      EqE2[i] = Eq_hold[1];
      EqE3[i] = Eq_hold[2];
      EqE4[i] = Eq_hold[3];
      EqW1[i] = Eq_hold[4];
      EqW2[i] = Eq_hold[5];
      EqW3[i] = Eq_hold[6];
      EqW4[i] = Eq_hold[7];
      //simFileIn >> EqE1[i] >> EqE2[i] >> EqE3[i] >> EqE4[i]
      //	>> EqW1[i] >> EqW2[i] >> EqW3[i] >> EqW4[i];
      //cout << EqE1[i] << " " << EqE2[i] << " " << EqE3[i] << " " << EqE4[i] << " "
      //   << EqW1[i] << " " << EqW2[i] << " " << EqW3[i] << " " << EqW4[i] << endl;
    }
    else {
      cout << "The simulation and data source peak files don't match. You need to re-make the source calibration files for both sim and data\n";
      exit(0);
    }

    if (filein.fail()) break;
    i++;
  }
  Int_t num = i;
  cout << "Number of data points: " << num << endl;

  /* // Load the smeared EQ values which are different for each PMT and source
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
    }*/

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

  // Filling new vectors for each PMT with ADC values only when the PMT is 
  // usable
  vector<int> runE1,runE2,runE3,runE4,runW1,runW2,runW3,runW4;
  vector<Double_t> EQE1,EQE2,EQE3,EQE4,EQW1,EQW2,EQW3,EQW4;
  vector<string> nameE1,nameE2,nameE3,nameE4,nameW1,nameW2,nameW3,nameW4;
  vector<Double_t> ADCE1, ADCE2, ADCE3, ADCE4;
  vector<Double_t> ADCW1, ADCW2, ADCW3, ADCW4;
  vector<Double_t> ResE1, ResE2, ResE3, ResE4;
  vector<Double_t> ResW1, ResW2, ResW3, ResW4;
  // Fill run, adc, and eQ vectors
  for (Int_t i=0; i<num;i++) {
    UInt_t runPos = find_vec_location_int(pmtRun,run[i]);
    cout << "Found run " <<  run[i] << " in PMT Quality\n";
    /*int src_hold;
    if (sourceName[i]=="Ce") src_hold = 0;
    else if (sourceName[i]=="Sn") src_hold=1;
    else if (sourceName[i]=="Bi2") src_hold=2;
    else if (sourceName[i]=="Bi1") src_hold=3;*/

    if (pmtQuality[runPos][0]) {
      runE1.push_back(run[i]);
      EQE1.push_back(EqE1[i]);
      nameE1.push_back(sourceName[i]);
      ADCE1.push_back(adcE1[i]);
    }
    if (pmtQuality[runPos][1]) {
      runE2.push_back(run[i]);
      EQE2.push_back(EqE2[i]);
      nameE2.push_back(sourceName[i]);
      ADCE2.push_back(adcE2[i]);
    }
    if (pmtQuality[runPos][2]) {
      runE3.push_back(run[i]);
      EQE3.push_back(EqE3[i]);
      nameE3.push_back(sourceName[i]);
      ADCE3.push_back(adcE3[i]);
    }
    if (pmtQuality[runPos][3]) {
      runE4.push_back(run[i]);
      EQE4.push_back(EqE4[i]);
      nameE4.push_back(sourceName[i]);
      ADCE4.push_back(adcE4[i]);
    }
    if (pmtQuality[runPos][4]) {
      runW1.push_back(run[i]);
      EQW1.push_back(EqW1[i]);
      nameW1.push_back(sourceName[i]);
      ADCW1.push_back(adcW1[i]);
    }
    if (pmtQuality[runPos][5]) {
      runW2.push_back(run[i]);
      EQW2.push_back(EqW2[i]);
      nameW2.push_back(sourceName[i]);
      ADCW2.push_back(adcW2[i]);
    }
    if (pmtQuality[runPos][6]) {
      runW3.push_back(run[i]);
      EQW3.push_back(EqW3[i]);
      nameW3.push_back(sourceName[i]);
      ADCW3.push_back(adcW3[i]);
    }
    if (pmtQuality[runPos][7]) {
      runW4.push_back(run[i]);
      EQW4.push_back(EqW4[i]);
      nameW4.push_back(sourceName[i]);
      ADCW4.push_back(adcW4[i]);
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
    

  //File to hold the linearity curves for each calibration run period
  sprintf(temp,"../linearity_curves/lin_curves_srcCal_Period_%i.dat",calibrationPeriod);
  ofstream linCurves(temp);

  // Fit function
  TF1 *fitADC = new TF1("fitADC", "([0] + [1]*x + [2]*x*x)*(0.5+0.5*TMath::TanH((x-[4])/[5]))+([3]*x)*(0.5-0.5*TMath::TanH((x-[4])/[5]))", 0., 2500.0);
  fitADC->SetParameter(0, 0.0);
  fitADC->SetParameter(1, 1.0);
  fitADC->SetParameter(2, 0.0);
  fitADC->SetParameter(3, 1.0);
  fitADC->SetParameter(4, 50.);
  fitADC->SetParameter(5, 1.);
  fitADC->FixParameter(5, 75.0);
  //fitADC->SetParLimits(5, 1., 20.);
  
  fitADC->SetParLimits(2, -0.0004, 0.0004);
  fitADC->SetParLimits(3, 0.1, 2.);

  fitADC->SetNpx(100000);
  fitADC->SetLineColor(2);

  Double_t offset=0., slope=0., quad=0., cubic=0.; //Fit Parameters
  Double_t slopeToOrigin = 0.; //This is the slope calculated from the low end of the fit range to the origin
  Double_t lowFitThreshold = 0.; //This is the low end of the quadratic fit region defined to be 
                                 // (1/2)*(Average Ce Peak in ADC).
  Double_t highFitThreshold = 0.; // This is the upper end of the quadratic fit
  Double_t lowEnSlope = 0., shiftOfTanH = 0., spreadOfTanH = 0.;

  // East 1

  sprintf(temp,"../residuals/residuals_East_runPeriod_%i_PMTE1.dat",calibrationPeriod);
  ofstream oFileE1(temp);
  vector <Double_t> fitEQ_E1(runE1.size(),0);

  if (runE1.size()>0 && (std::find(nameE1.begin(),nameE1.end(),"Ce")!=nameE1.end() && std::find(nameE1.begin(),nameE1.end(),"Sn")!=nameE1.end() && std::find(nameE1.begin(),nameE1.end(),"Bi1")!=nameE1.end())) {
 
    //Calculating the average of the lower Ce ADC values and assigning a value to lowFitThreshold
    Double_t sumADC_Ce = 0., sumADC_Bi = 0., sumEQ_Ce=0.; 
    Int_t entries1 = 0, entries2 = 0;;
    for (UInt_t ii=0; ii<runE1.size(); ii++) {
      if (nameE1[ii]=="Ce") {
	sumADC_Ce+=ADCE1[ii];
	sumEQ_Ce+=EQE1[ii];
	entries1++;
      }
      if (nameE1[ii]=="Bi1") {
	sumADC_Bi+=ADCE1[ii];
	entries2++;
      }
    }
    
    //Calculating the slope of a line from the origin to the mean of the Ce peaks
    Double_t linSlope = (sumEQ_Ce/sumADC_Ce);
    fitADC->FixParameter(3,linSlope);
    fitADC->FixParameter(4,sumADC_Ce/entries1);
    //lowFitThreshold = sumADC_Ce/entries1;
    highFitThreshold = 1.5*sumADC_Bi/entries2;
    //cout << lowFitThreshold << endl;
    //cout << highFitThreshold << endl;

    c1 = new TCanvas("c1", "c1");
    c1->SetLogy(0);

    TGraphErrors *grE1 = new TGraphErrors(runE1.size(),&ADCE1[0],&EQE1[0],err,err);
    grE1->SetTitle("");
    grE1->SetMarkerColor(1);
    grE1->SetLineColor(1);
    grE1->SetMarkerStyle(20);
    grE1->SetMarkerSize(0.75);
    grE1->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
    grE1->GetXaxis()->SetTitleOffset(1.2);
    grE1->GetXaxis()->CenterTitle();
    grE1->GetYaxis()->SetTitle("E_{Q} [keV]");
    grE1->GetYaxis()->SetTitleOffset(1.6);
    grE1->GetYaxis()->CenterTitle();
    grE1->GetXaxis()->SetLimits(0.0,2400.);
    grE1->SetMinimum(0.0);
    grE1->SetMaximum(1000.0);
    grE1->Draw("AP");

    grE1->Fit("fitADC", "","", 0., highFitThreshold);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);
    lowEnSlope = fitADC->GetParameter(3);
    shiftOfTanH = fitADC->GetParameter(4);
    spreadOfTanH = fitADC->GetParameter(5);
    slopeToOrigin = (offset + slope*lowFitThreshold + quad*lowFitThreshold*lowFitThreshold + cubic*lowFitThreshold*lowFitThreshold*lowFitThreshold)/lowFitThreshold;
    linCurves << offset << " " << slope << " " << quad << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;
  
    Double_t x1_text = 1200;
    Double_t y1_text = 100;

    TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("East PMT 1");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();

    // Calculate residuals in [keV]
    
    for (int j=0; j<runE1.size(); j++) {
    
      fitEQ_E1[j]    = fitADC->Eval(ADCE1[j]);//offset + slope*ADCE1[j] + quad*ADCE1[j]*ADCE1[j] + cubic*ADCE1[j]*ADCE1[j]*ADCE1[j];
      if (nameE1[j]=="Ce") {
	ResE1[j] = fitEQ_E1[j] - EQE1[j];
	oFileE1 << "Ce_East" << " " << runE1[j] << " " << ResE1[j] << endl;
	if (sqrt(ResE1[j]*ResE1[j])>0.025*EQE1[j]) cout << "Ce_East" << " PMT1 " << runE1[j] << " " << ResE1[j] << endl;
	//ResE1[j] = (fitEQ - peakCe)/peakCe * 100.;
      }
      else if (nameE1[j]=="Sn") {
	ResE1[j] = fitEQ_E1[j] - EQE1[j];
	//ResE1[j] = (fitEQ - peakSn)/peakSn * 100.;
	oFileE1 << "Sn_East" << " " << runE1[j] << " " << ResE1[j] << endl;
	if (sqrt(ResE1[j]*ResE1[j])>0.025*EQE1[j]) cout << "Sn_East" << " " << "PMT1 " << runE1[j] << " " << ResE1[j] << endl;
      }
      else if (nameE1[j]=="Bi1") {
	ResE1[j] = fitEQ_E1[j] - EQE1[j];
	//ResE1[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE1 << "Bi1_East" << " " << runE1[j] << " " << ResE1[j] << endl;
	if (sqrt(ResE1[j]*ResE1[j])>0.025*EQE1[j]) cout << "Bi1_East" << " " << "PMT1 " << runE1[j] << " " << ResE1[j] << endl;
 
      }
    }

    // East 1 Residuals
    c1r = new TCanvas("c1r", "c1r");
    c1r->SetLogy(0);

    TGraphErrors *grE1r = new TGraphErrors(runE1.size(),&ADCE1[0],&ResE1[0],err,err);
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
    grE1r->GetXaxis()->SetLimits(0.0,2400.);
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

  else { linCurves << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
    oFileE1 << "PMT NOT USABLE";}
  oFileE1.close();

  // East 2

  sprintf(temp,"../residuals/residuals_East_runPeriod_%i_PMTE2.dat",calibrationPeriod);
  ofstream oFileE2(temp);
  vector <Double_t> fitEQ_E2(runE2.size(),0);

  if (runE2.size()>0 && (std::find(nameE2.begin(),nameE2.end(),"Ce")!=nameE2.end() && std::find(nameE2.begin(),nameE2.end(),"Sn")!=nameE2.end() && std::find(nameE2.begin(),nameE2.end(),"Bi1")!=nameE2.end())) {

    //Calculating the average of the lower Ce ADC values and assigning a value to lowFitThreshold
    Double_t sumADC_Ce = 0., sumADC_Bi = 0., sumEQ_Ce=0.; 
    Int_t entries1 = 0, entries2 = 0;;
    for (UInt_t ii=0; ii<runE2.size(); ii++) {
      if (nameE2[ii]=="Ce") {
	sumADC_Ce+=ADCE2[ii];
	sumEQ_Ce+=EQE2[ii];
	entries1++;
      }
      if (nameE2[ii]=="Bi1") {
	sumADC_Bi+=ADCE2[ii];
	entries2++;
      }
    }
    
    //Calculating the slope of a line from the origin to the mean of the Ce peaks
    Double_t linSlope = (sumEQ_Ce/sumADC_Ce);
    fitADC->FixParameter(3,linSlope);
    fitADC->FixParameter(4,sumADC_Ce/entries1);
    //lowFitThreshold = sumADC_Ce/entries1;
    highFitThreshold = 1.5*sumADC_Bi/entries2;
    //cout << lowFitThreshold << endl;
    //cout << highFitThreshold << endl;
    

    c2 = new TCanvas("c2", "c2");
    c2->SetLogy(0);

    TGraphErrors *grE2 = new TGraphErrors(runE2.size(),&ADCE2[0],&EQE2[0],err,err);
    grE2->SetTitle("");
    grE2->SetMarkerColor(1);
    grE2->SetLineColor(1);
    grE2->SetMarkerStyle(20);
    grE2->SetMarkerSize(0.75);
    grE2->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
    grE2->GetXaxis()->SetTitleOffset(1.2);
    grE2->GetXaxis()->CenterTitle();
    grE2->GetYaxis()->SetTitle("E_{Q} [keV]");
    grE2->GetYaxis()->SetTitleOffset(1.6);
    grE2->GetYaxis()->CenterTitle();
    grE2->GetXaxis()->SetLimits(0.0,2400.);
    grE2->SetMinimum(0.0);
    grE2->SetMaximum(1000.0);
    grE2->Draw("AP");

    grE2->Fit("fitADC", "Q","",0.,highFitThreshold);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);
    lowEnSlope = fitADC->GetParameter(3);
    shiftOfTanH = fitADC->GetParameter(4);
    spreadOfTanH = fitADC->GetParameter(5);
    slopeToOrigin = (offset + slope*lowFitThreshold + quad*lowFitThreshold*lowFitThreshold + cubic*lowFitThreshold*lowFitThreshold*lowFitThreshold)/lowFitThreshold;
    linCurves << offset << " " << slope << " " << quad << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;

    Double_t x1_text = 1200;
    Double_t y1_text = 100;

    TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("East PMT 2");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();

    // Calculate residuals in [keV]
    
    for (int j=0; j<runE2.size(); j++) {
      fitEQ_E2[j]    = fitADC->Eval(ADCE2[j]);//offset + slope*ADCE2[j] + quad*ADCE2[j]*ADCE2[j] + cubic*ADCE2[j]*ADCE2[j]*ADCE2[j];
      if (nameE2[j]=="Ce") {
	ResE2[j] = fitEQ_E2[j] - EQE2[j];
	//ResE2[j] = (fitEQ - peakCe)/peakCe * 100.;
	oFileE2 << "Ce_East" << " " << runE2[j] << " " << ResE2[j] << endl;
	if (sqrt(ResE2[j]*ResE2[j])>0.025*EQE2[j]) cout<< "Ce_East" << " " << "PMT2 " << runE2[j] << " " << ResE2[j] << endl;  
      }
      else if (nameE2[j]=="Sn") {
	ResE2[j] = fitEQ_E2[j] - EQE2[j];
	//ResE2[j] = (fitEQ - peakSn)/peakSn * 100.;
	oFileE2 << "Sn_East" << " " << runE2[j] << " " << ResE2[j] << endl;
	if (sqrt(ResE2[j]*ResE2[j])>0.025*EQE2[j]) cout<< "Sn_East" << " " << "PMT2 " << runE2[j] << " " << ResE2[j] << endl;  
      }
      else if (nameE2[j]=="Bi1") {
	ResE2[j] = fitEQ_E2[j] - EQE2[j];
	//ResE2[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE2 << "Bi1_East" << " " << runE2[j] << " " << ResE2[j] << endl;
	if (sqrt(ResE2[j]*ResE2[j])>0.025*EQE2[j]) cout<< "Bi1_East" << " " << "PMT2 " << runE2[j] << " " << ResE2[j] << endl;  
      }
    }

    // East 2 residuals
    c2r = new TCanvas("c2r", "c2r");
    c2r->SetLogy(0);

    TGraphErrors *grE2r = new TGraphErrors(runE2.size(),&ADCE2[0],&ResE2[0],err,err);
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
    grE2r->GetXaxis()->SetLimits(0.0,2400.);
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

  else { linCurves << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
    oFileE2 << "PMT NOT USABLE";}
  oFileE2.close();

  // East 3

  sprintf(temp,"../residuals/residuals_East_runPeriod_%i_PMTE3.dat",calibrationPeriod);
  ofstream oFileE3(temp);
  vector <Double_t> fitEQ_E3(runE3.size(),0);

  if (runE3.size()>0 && (std::find(nameE3.begin(),nameE3.end(),"Ce")!=nameE3.end() && std::find(nameE3.begin(),nameE3.end(),"Sn")!=nameE3.end() && std::find(nameE3.begin(),nameE3.end(),"Bi1")!=nameE3.end())) {

    //Calculating the average of the lower Ce ADC values and assigning a value to lowFitThreshold
    Double_t sumADC_Ce = 0., sumADC_Bi = 0., sumEQ_Ce=0.; 
    Int_t entries1 = 0, entries2 = 0;;
    for (UInt_t ii=0; ii<runE3.size(); ii++) {
      if (nameE3[ii]=="Ce") {
	sumADC_Ce+=ADCE3[ii];
	sumEQ_Ce+=EQE3[ii];
	entries1++;
      }
      if (nameE3[ii]=="Bi1") {
	sumADC_Bi+=ADCE3[ii];
	entries2++;
      }
    }
    
    //Calculating the slope of a line from the origin to the mean of the Ce peaks
    Double_t linSlope = (sumEQ_Ce/sumADC_Ce);
    fitADC->FixParameter(3,linSlope);
    fitADC->FixParameter(4,sumADC_Ce/entries1);
    //lowFitThreshold = sumADC_Ce/entries1;
    highFitThreshold = 1.5*sumADC_Bi/entries2;
    //cout << lowFitThreshold << endl;
    //cout << highFitThreshold << endl;

    c3 = new TCanvas("c3", "c3");
    c3->SetLogy(0);

    TGraphErrors *grE3 = new TGraphErrors(runE3.size(),&ADCE3[0],&EQE3[0],err,err);
    grE3->SetTitle("");
    grE3->SetMarkerColor(1);
    grE3->SetLineColor(1);
    grE3->SetMarkerStyle(20);
    grE3->SetMarkerSize(0.75);
    grE3->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
    grE3->GetXaxis()->SetTitleOffset(1.2);
    grE3->GetXaxis()->CenterTitle();
    grE3->GetYaxis()->SetTitle("E_{Q} [keV]");
    grE3->GetYaxis()->SetTitleOffset(1.6);
    grE3->GetYaxis()->CenterTitle();
    grE3->GetXaxis()->SetLimits(0.0,2000.0);
    grE3->SetMinimum(0.0);
    grE3->SetMaximum(1000.0);
    grE3->Draw("AP");

    grE3->Fit("fitADC", "","",0.,highFitThreshold);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);
    lowEnSlope = fitADC->GetParameter(3);
    shiftOfTanH = fitADC->GetParameter(4);
    spreadOfTanH = fitADC->GetParameter(5);
    slopeToOrigin = (offset + slope*lowFitThreshold + quad*lowFitThreshold*lowFitThreshold + cubic*lowFitThreshold*lowFitThreshold*lowFitThreshold)/lowFitThreshold;
    linCurves << offset << " " << slope << " " << quad << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;

    Double_t x1_text = 1200;
    Double_t y1_text = 100;

    TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("East PMT 3");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();

    // Calculate residuals in [keV]
    
    for (int j=0; j<runE3.size(); j++) {
      fitEQ_E3[j]    = fitADC->Eval(ADCE3[j]);//offset + slope*ADCE3[j] + quad*ADCE3[j]*ADCE3[j] + cubic*ADCE3[j]*ADCE3[j]*ADCE3[j];
      if (nameE3[j]=="Ce") {
	ResE3[j] = fitEQ_E3[j] - EQE3[j];
	//ResE3[j] = (fitEQ - peakCe)/peakCe * 100.;
	oFileE3 << "Ce_East" << " " << runE3[j] << " " << ResE3[j] << endl;
	if (sqrt(ResE3[j]*ResE3[j])>0.025*EQE3[j]) cout<< "Ce_East" << " " << "PMT3 " << runE3[j] << " " << ResE3[j] << endl;  
      }
      else if (nameE3[j]=="Sn") {
	ResE3[j] = fitEQ_E3[j] - EQE3[j];
	//ResE3[j] = (fitEQ - peakSn)/peakSn * 100.;
	oFileE3 << "Sn_East" << " " << runE3[j] << " " << ResE3[j] << endl;
	if (sqrt(ResE3[j]*ResE3[j])>0.025*EQE3[j]) cout << "Sn_East" << " " << "PMT3 " << runE3[j] << " " << ResE3[j] << endl;  
      }
      else if (nameE3[j]=="Bi1") {
	ResE3[j] = fitEQ_E3[j] - EQE3[j];
	//ResE3[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE3 << "Bi1_East" << " " << runE3[j] << " " << ResE3[j] << endl;
	if (sqrt(ResE3[j]*ResE3[j])>0.025*EQE3[j]) cout << "Bi1_East" << " " << "PMT3 " << runE3[j] << " " << ResE3[j] << endl;  
      }
    }

    // East 3 residuals
    c3r = new TCanvas("c3r", "c3r");
    c3r->SetLogy(0);

    TGraphErrors *grE3r = new TGraphErrors(runE3.size(),&ADCE3[0],&ResE3[0],err,err);
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
    grE3r->GetXaxis()->SetLimits(0.0,2400.);
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

  else { linCurves << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
    oFileE3 << "PMT NOT USABLE";}
  oFileE3.close();

  // East 4

  sprintf(temp,"../residuals/residuals_East_runPeriod_%i_PMTE4.dat",calibrationPeriod);
  ofstream oFileE4(temp);
  vector <Double_t> fitEQ_E4(runE4.size(),0);

  if (runE4.size()>0 && (std::find(nameE4.begin(),nameE4.end(),"Ce")!=nameE4.end() && std::find(nameE4.begin(),nameE4.end(),"Sn")!=nameE4.end() && std::find(nameE4.begin(),nameE4.end(),"Bi1")!=nameE4.end())) {

    //Calculating the average of the lower Ce ADC values and assigning a value to lowFitThreshold
    Double_t sumADC_Ce = 0., sumADC_Bi = 0., sumEQ_Ce=0.; 
    Int_t entries1 = 0, entries2 = 0;;
    for (UInt_t ii=0; ii<runE4.size(); ii++) {
      if (nameE4[ii]=="Ce") {
	sumADC_Ce+=ADCE4[ii];
	sumEQ_Ce+=EQE4[ii];
	entries1++;
      }
      if (nameE4[ii]=="Bi1") {
	sumADC_Bi+=ADCE4[ii];
	entries2++;
      }
    }
    
    //Calculating the slope of a line from the origin to the mean of the Ce peaks
    Double_t linSlope = (sumEQ_Ce/sumADC_Ce);
    fitADC->FixParameter(3,linSlope);
    fitADC->FixParameter(4,sumADC_Ce/entries1);
    //lowFitThreshold = sumADC_Ce/entries1;
    highFitThreshold = 1.5*sumADC_Bi/entries2;
    //cout << lowFitThreshold << endl;
    //cout << highFitThreshold << endl;

    c4 = new TCanvas("c4", "c4");
    c4->SetLogy(0);

    TGraphErrors *grE4 = new TGraphErrors(runE4.size(),&ADCE4[0],&EQE4[0],err,err);
    grE4->SetTitle("");
    grE4->SetMarkerColor(1);
    grE4->SetLineColor(1);
    grE4->SetMarkerStyle(20);
    grE4->SetMarkerSize(0.75);
    grE4->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
    grE4->GetXaxis()->SetTitleOffset(1.2);
    grE4->GetXaxis()->CenterTitle();
    grE4->GetYaxis()->SetTitle("E_{Q} [keV]");
    grE4->GetYaxis()->SetTitleOffset(1.6);
    grE4->GetYaxis()->CenterTitle();
    grE4->GetXaxis()->SetLimits(0.0,2400.);
    grE4->SetMinimum(0.0);
    grE4->SetMaximum(1000.0);
    grE4->Draw("AP");

    grE4->Fit("fitADC", "","",0.,highFitThreshold);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);
    lowEnSlope = fitADC->GetParameter(3);
    shiftOfTanH = fitADC->GetParameter(4);
    spreadOfTanH = fitADC->GetParameter(5);
    slopeToOrigin = (offset + slope*lowFitThreshold + quad*lowFitThreshold*lowFitThreshold + cubic*lowFitThreshold*lowFitThreshold*lowFitThreshold)/lowFitThreshold;
    linCurves << offset << " " << slope << " " << quad << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;
    
    Double_t x1_text = 1200;
    Double_t y1_text = 100;

    TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("East PMT 4");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();

    // Calculate residuals in [keV]
    
    for (int j=0; j<runE4.size(); j++) {
      fitEQ_E4[j]    = fitADC->Eval(ADCE4[j]);//offset + slope*ADCE4[j] + quad*ADCE4[j]*ADCE4[j] + cubic*ADCE4[j]*ADCE4[j]*ADCE4[j];
      if (nameE4[j]=="Ce") {
	ResE4[j] = fitEQ_E4[j] - EQE4[j];
	//ResE4[j] = (fitEQ - peakCe)/peakCe * 100.;
	oFileE4 << "Ce_East" << " " << runE4[j] << " " << ResE4[j] << endl;
	if (sqrt(ResE4[j]*ResE4[j])>0.025*EQE4[j]) cout << "Ce_East" << " " << "PMT4 " << runE4[j] << " " << ResE4[j] << endl;  
      }
      else if (nameE4[j]=="Sn") {
	ResE4[j] = fitEQ_E4[j] - EQE4[j];
	//ResE4[j] = (fitEQ - peakSn)/peakSn * 100.;
	oFileE4 << "Sn_East" << " " << runE4[j] << " " << ResE4[j] << endl;
	if (sqrt(ResE4[j]*ResE4[j])>0.025*EQE4[j]) cout<< "Sn_East" << " " << "PMT4 " << runE4[j] << " " << ResE4[j] << endl;  
      }
      else if (nameE4[j]=="Bi1") {
	ResE4[j] = fitEQ_E4[j] - EQE4[j];
	//ResE4[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE4 << "Bi1_East" << " " << runE4[j] << " " << ResE4[j] << endl;
	if (sqrt(ResE4[j]*ResE4[j])>0.025*EQE4[j]) cout<< "Bi1_East" << " " << "PMT4 " << runE4[j] << " " << ResE4[j] << endl;  
      }
    }

    // East 4 residuals
    c4r = new TCanvas("c4r", "c4r");
    c4r->SetLogy(0);

    TGraphErrors *grE4r = new TGraphErrors(runE4.size(),&ADCE4[0],&ResE4[0],err,err);
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
    grE4r->GetXaxis()->SetLimits(0.0,2400.);
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

  else { linCurves << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
    oFileE4 << "PMT NOT USABLE";}
  oFileE4.close();

  ///////////////////////////////////////////////////////////////////////
  // Calcuting the weighted mean of the total energy of the East side  //
  ///////////////////////////////////////////////////////////////////////

  sprintf(temp,"../residuals/residuals_East_runPeriod_%i.dat",calibrationPeriod);
  ofstream oFileE(temp);

  cout << "CALCULATING EAST RESIDUALS" << endl;
  vector <Double_t> Etrue_East(num,0);
  vector <Double_t> x_East(num,0);
  vector <Double_t> res_East(num,0);

  std::vector<Int_t>::iterator E1 = runE1.begin();
  std::vector<Double_t>::iterator E1_EQ = EQE1.begin();
  std::vector<Double_t>::iterator E1_fitEQ = fitEQ_E1.begin();
  std::vector<Int_t>::iterator E2 = runE2.begin();
  std::vector<Double_t>::iterator E2_EQ = EQE2.begin();
  std::vector<Double_t>::iterator E2_fitEQ = fitEQ_E2.begin();
  std::vector<Int_t>::iterator E3 = runE3.begin();
  std::vector<Double_t>::iterator E3_EQ = EQE3.begin();
  std::vector<Double_t>::iterator E3_fitEQ = fitEQ_E3.begin();
  std::vector<Int_t>::iterator E4 = runE4.begin();
  std::vector<Double_t>::iterator E4_EQ = EQE4.begin();
  std::vector<Double_t>::iterator E4_fitEQ = fitEQ_E4.begin();
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
    UInt_t pmtVecLocation = find_vec_location_int(pmtRun,run[j]);

    if (pmtQuality[pmtVecLocation][0] && *E1==run[j] && *nmE1==sourceName[j]) {
      Double_t N = adcE1[j]*nPE_per_channel[pmtVecLocation][0];
      Double_t f = sqrt(N)/N;
      ResidualE1 = *RE1;
      //Energy1 = EQ2Etrue[0][0]+EQ2Etrue[0][1]*(*E1_fitEQ)+EQ2Etrue[0][2]*(*E1_fitEQ)*(*E1_fitEQ);
      Energy1 = E1_fitEQ;
      weight1 = 1/(Energy1*Energy1*f*f);
      std::advance(E1,1);
      std::advance(E1_EQ,1);
      std::advance(E1_fitEQ,1);
      std::advance(nmE1,1);
      std::advance(RE1,1);
    }
    else {weight1=0.; Energy1=0.; ResidualE1=0.;}

    if (pmtQuality[pmtVecLocation][1] && *E2==run[j] && *nmE2==sourceName[j]) {
      Double_t N = adcE2[j]*nPE_per_channel[pmtVecLocation][1];
      Double_t f = sqrt(N)/N;
      ResidualE2 = *RE2;
      //Energy2 = EQ2Etrue[1][0]+EQ2Etrue[1][1]*(*E2_fitEQ)+EQ2Etrue[1][2]*(*E2_fitEQ)*(*E2_fitEQ);
      Energy2 = E2_fitEQ;
      weight2 = 1/(Energy2*Energy2*f*f);
      std::advance(E2,1);
      std::advance(E2_EQ,1);
      std::advance(E2_fitEQ,1);
      std::advance(nmE2,1);
      std::advance(RE2,1);
    }
    else {weight2=0.; Energy2=0.; ResidualE2=0.;}

    if (pmtQuality[pmtVecLocation][2] && *E3==run[j] && *nmE3==sourceName[j]) {
      Double_t N = adcE3[j]*nPE_per_channel[pmtVecLocation][2];
      Double_t f = sqrt(N)/N;
      ResidualE4 = *RE3;
      //Energy3 = EQ2Etrue[2][0]+EQ2Etrue[2][1]*(*E3_fitEQ)+EQ2Etrue[2][2]*(*E3_fitEQ)*(*E3_fitEQ);;
      Energy3 = E3_fitEQ;
      weight3 = 1/(Energy3*Energy3*f*f);
      std::advance(E3,1);
      std::advance(E3_EQ,1);
      std::advance(E3_fitEQ,1);
      std::advance(nmE3,1);
      std::advance(RE3,1);
    }
    else {weight3=0.; Energy3=0.; ResidualE3=0.;}

    if (pmtQuality[pmtVecLocation][3] && *E4==run[j] && *nmE4==sourceName[j]) {
      Double_t N = adcE4[j]*nPE_per_channel[pmtVecLocation][3];
      Double_t f = sqrt(N)/N;
      ResidualE4 = *RE4;
      //Energy4 = EQ2Etrue[3][0]+EQ2Etrue[3][1]*(*E4_fitEQ)+EQ2Etrue[3][2]*(*E4_fitEQ)*(*E4_fitEQ);
      Energy4 = E4_fitEQ;
      weight4 = 1/(Energy4*Energy4*f*f);
      std::advance(E4,1);
      std::advance(E4_EQ,1);
      std::advance(E4_fitEQ,1);
      std::advance(nmE4,1);
      std::advance(RE4,1);
    }
    else {weight4=0.; Energy4=0.; ResidualE4=0.;}
  
    Etrue_East[j] = (weight1*Energy1+weight2*Energy2+weight3*Energy3+weight4*Energy4)/(weight1+weight2+weight3+weight4);
    res_East[j] = (weight1*ResidualE1+weight2*ResidualE2+weight3*ResidualE3+weight4*ResidualE4)/(weight1+weight2+weight3+weight4);
    if (sourceName[j]=="Ce") {
      //res_East[j] = Etrue_East[j] - peakCe;
      x_East[j] = peakCe_EQ;
      oFileE << "Ce_East" << " " <<  run[j] << " " << res_East[j] << endl;
    }
    else if (sourceName[j]=="Sn") {
      //res_East[j] = Etrue_East[j] - peakSn;
      x_East[j] = peakSn_EQ;
      oFileE << "Sn_East" << " " <<  run[j] << " " << res_East[j] << endl;
      /*cout << Energy1 << " " << weight1 << " "
	   << Energy2 << " " << weight2 << " "
	   << Energy3 << " " << weight3 << " " 
	   << Energy4 << " " << weight4 << endl;*/
    }
    else if (sourceName[j]=="Bi2") {
      //res_East[j] = Etrue_East[j] - peakBiLow;
      x_East[j] = peakBiLow_EQ;
    }
    else if (sourceName[j]=="Bi1") {
      //res_East[j] = Etrue_East[j] - peakBiHigh;
      x_East[j] = peakBiHigh_EQ;
      oFileE << "Bi1_East" << " " <<  run[j] << " " << res_East[j] << endl;
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

  sprintf(temp,"../residuals/residuals_West_runPeriod_%i_PMTW1.dat",calibrationPeriod);
  ofstream oFileW1(temp);
  vector <Double_t> fitEQ_W1(runW1.size(),0);

  if (runW1.size()>0 && (std::find(nameW1.begin(),nameW1.end(),"Ce")!=nameW1.end() && std::find(nameW1.begin(),nameW1.end(),"Sn")!=nameW1.end() && std::find(nameW1.begin(),nameW1.end(),"Bi1")!=nameW1.end())) {

    //Calculating the average of the lower Ce ADC values and assigning a value to lowFitThreshold
    Double_t sumADC_Ce = 0., sumADC_Bi = 0., sumEQ_Ce=0.; 
    Int_t entries1 = 0, entries2 = 0;;
    for (UInt_t ii=0; ii<runW1.size(); ii++) {
      if (nameW1[ii]=="Ce") {
	sumADC_Ce+=ADCW1[ii];
	sumEQ_Ce+=EQW1[ii];
	entries1++;
      }
      if (nameW1[ii]=="Bi1") {
	sumADC_Bi+=ADCW1[ii];
	entries2++;
      }
    }
    
    //Calculating the slope of a line from the origin to the mean of the Ce peaks
    Double_t linSlope = (sumEQ_Ce/sumADC_Ce);
    fitADC->FixParameter(3,linSlope);
    fitADC->FixParameter(4,sumADC_Ce/entries1);
    //lowFitThreshold = sumADC_Ce/entries1;
    highFitThreshold = 1.5*sumADC_Bi/entries2;
    //cout << lowFitThreshold << endl;
    //cout << highFitThreshold << endl;

    cW1 = new TCanvas("cW1", "cW1");
    cW1->SetLogy(0);

    TGraphErrors *grW1 = new TGraphErrors(runW1.size(),&ADCW1[0],&EQW1[0],err,err);
    grW1->SetTitle("");
    grW1->SetMarkerColor(1);
    grW1->SetLineColor(1);
    grW1->SetMarkerStyle(20);
    grW1->SetMarkerSize(0.75);
    grW1->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
    grW1->GetXaxis()->SetTitleOffset(1.2);
    grW1->GetXaxis()->CenterTitle();
    grW1->GetYaxis()->SetTitle("E_{Q} [keV]");
    grW1->GetYaxis()->SetTitleOffset(1.6);
    grW1->GetYaxis()->CenterTitle();
    grW1->GetXaxis()->SetLimits(0.0,2400.);
    grW1->SetMinimum(0.0);
    grW1->SetMaximum(1000.0);
    grW1->Draw("AP");

    grW1->Fit("fitADC", "","",0.,highFitThreshold);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);
    lowEnSlope = fitADC->GetParameter(3);
    shiftOfTanH = fitADC->GetParameter(4);
    spreadOfTanH = fitADC->GetParameter(5);
    slopeToOrigin = (offset + slope*lowFitThreshold + quad*lowFitThreshold*lowFitThreshold + cubic*lowFitThreshold*lowFitThreshold*lowFitThreshold)/lowFitThreshold;
    linCurves << offset << " " << slope << " " << quad << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;

    Double_t x1_text = 1200;
    Double_t y1_text = 100;

    TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("West PMT 1");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();
  
    // Calculate residuals in [keV]
    //Double_t fitEQ_W1[num];
    for (int j=0; j<runW1.size(); j++) {
      fitEQ_W1[j]    = fitADC->Eval(ADCW1[j]);//offset + slope*ADCW1[j] + quad*ADCW1[j]*ADCW1[j] + cubic*ADCW1[j]*ADCW1[j]*ADCW1[j];
      if (nameW1[j]=="Ce") {
	ResW1[j] = fitEQ_W1[j] - EQW1[j];
	oFileW1 << "Ce_West" << " " << runW1[j] << " " << ResW1[j] << endl;
	if (sqrt(ResW1[j]*ResW1[j])>0.025*EQW1[j]) cout<< "Ce_West" << " " << "PMT1 " << runW1[j] << " " << ResW1[j] << endl;  
      }
      else if (nameW1[j]=="Sn") {
	ResW1[j] = fitEQ_W1[j] - EQW1[j];
	oFileW1 << "Sn_West" << " " <<  runW1[j] << " " << ResW1[j] << endl;
	if (sqrt(ResW1[j]*ResW1[j])>0.025*EQW1[j]) cout<< "Sn_West" << " " << "PMT1 " << runW1[j] << " " << ResW1[j] << endl;  
      }
      else if (nameW1[j]=="Bi1") {
	ResW1[j] = fitEQ_W1[j] - EQW1[j];
	oFileW1 << "Bi1_West" << " " << runW1[j] << " " << ResW1[j] << endl;
	if (sqrt(ResW1[j]*ResW1[j])>0.025*EQW1[j]) cout<< "Bi1_West" << " " << "PMT1 " << runW1[j] << " " << ResW1[j] << endl;  
      }
    }

    // West 1 residuals
    cW1r = new TCanvas("cW1r", "cW1r");
    cW1r->SetLogy(0);

    TGraphErrors *grW1r = new TGraphErrors(runW1.size(),&ADCW1[0],&ResW1[0],err,err);
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
    grW1r->GetXaxis()->SetLimits(0.0,2400.);
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

  else { linCurves << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
    oFileW1 << "PMT NOT USABLE";}
  oFileW1.close();

  // West 2

  sprintf(temp,"../residuals/residuals_West_runPeriod_%i_PMTW2.dat",calibrationPeriod);
  ofstream oFileW2(temp);
  vector <Double_t> fitEQ_W2(runW2.size(),0);

  if (runW2.size()>0 && (std::find(nameW2.begin(),nameW2.end(),"Ce")!=nameW2.end() && std::find(nameW2.begin(),nameW2.end(),"Sn")!=nameW2.end() && std::find(nameW2.begin(),nameW2.end(),"Bi1")!=nameW2.end())) {
    //Calculating the average of the lower Ce ADC values and assigning a value to lowFitThreshold
    Double_t sumADC_Ce = 0., sumADC_Bi = 0., sumEQ_Ce=0.; 
    Int_t entries1 = 0, entries2 = 0;;
    for (UInt_t ii=0; ii<runW2.size(); ii++) {
      if (nameW2[ii]=="Ce") {
	sumADC_Ce+=ADCW2[ii];
	sumEQ_Ce+=EQW2[ii];
	entries1++;
      }
      if (nameW2[ii]=="Bi1") {
	sumADC_Bi+=ADCW2[ii];
	entries2++;
      }
    }

    //Calculating the slope of a line from the origin to the mean of the Ce peaks
    Double_t linSlope = (sumEQ_Ce/sumADC_Ce);
    fitADC->FixParameter(3,linSlope);
    fitADC->FixParameter(4,sumADC_Ce/entries1);
    //lowFitThreshold = sumADC_Ce/entries1;
    highFitThreshold = 1.5*sumADC_Bi/entries2;
    //cout << lowFitThreshold << endl;
    //cout << highFitThreshold << endl;


    cW2 = new TCanvas("cW2", "cW2");
    cW2->SetLogy(0);

    TGraphErrors *grW2 = new TGraphErrors(runW2.size(),&ADCW2[0],&EQW2[0],err,err);
    grW2->SetTitle("");
    grW2->SetMarkerColor(1);
    grW2->SetLineColor(1);
    grW2->SetMarkerStyle(20);
    grW2->SetMarkerSize(0.75);
    grW2->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
    grW2->GetXaxis()->SetTitleOffset(1.2);
    grW2->GetXaxis()->CenterTitle();
    grW2->GetYaxis()->SetTitle("E_{Q} [keV]");
    grW2->GetYaxis()->SetTitleOffset(1.6);
    grW2->GetYaxis()->CenterTitle();
    grW2->GetXaxis()->SetLimits(0.0,2400.);
    grW2->SetMinimum(0.0);
    grW2->SetMaximum(1000.0);
    grW2->Draw("AP");

    grW2->Fit("fitADC", "","",0.,highFitThreshold);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);
    lowEnSlope = fitADC->GetParameter(3);
    shiftOfTanH = fitADC->GetParameter(4);
    spreadOfTanH = fitADC->GetParameter(5);
    slopeToOrigin = (offset + slope*lowFitThreshold + quad*lowFitThreshold*lowFitThreshold + cubic*lowFitThreshold*lowFitThreshold*lowFitThreshold)/lowFitThreshold;
    linCurves << offset << " " << slope << " " << quad << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;

    Double_t x1_text = 1200;
    Double_t y1_text = 100;

    TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("West PMT 2");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();

    // Calculate residuals in [keV]
    //Double_t fitEQ_W2[num];
    for (int j=0; j<runW2.size(); j++) {
      fitEQ_W2[j]    = fitADC->Eval(ADCW2[j]);//offset + slope*ADCW2[j] + quad*ADCW2[j]*ADCW2[j] + cubic*ADCW2[j]*ADCW2[j]*ADCW2[j];
      if (nameW2[j]=="Ce") {
	ResW2[j] = fitEQ_W2[j] - EQW2[j];
	oFileW2 << "Ce_West" << " " << runW2[j] << " " << ResW2[j] << endl;
	if (sqrt(ResW2[j]*ResW2[j])>0.025*EQW2[j]) cout<< "Ce_West" << " " << "PMT2 " << runW2[j] << " " << ResW2[j] << endl;  
      }
      else if (nameW2[j]=="Sn") {
	ResW2[j] = fitEQ_W2[j] - EQW2[j];
	oFileW2 << "Sn_West" << " " << runW2[j] << " " << ResW2[j] << endl;
	if (sqrt(ResW2[j]*ResW2[j])>0.025*EQW2[j]) cout<< "Sn_West" << " " << "PMT2 " << runW2[j] << " " << ResW2[j] << endl;  
      }
      else if (nameW2[j]=="Bi1") {
	ResW2[j] = fitEQ_W2[j] - EQW2[j];
	oFileW2 << "Bi1_West" << " " << runW2[j] << " " << ResW2[j] << endl;
	if (sqrt(ResW2[j]*ResW2[j])>0.025*EQW2[j]) cout<< "Bi1_West" << " " << "PMT2 " << runW2[j] << " " << ResW2[j] << endl;  
      }
    }

    // West 2 residuals
    cW2r = new TCanvas("cW2r", "cW2r");
    cW2r->SetLogy(0);

    TGraphErrors *grW2r = new TGraphErrors(runW2.size(),&ADCW2[0],&ResW2[0],err,err);
    grW2r->SetTitle("");
    grW2r->SetMarkerColor(1);
    grW2r->SetLineColor(1);
    grW2r->SetMarkerStyle(20);
    grW2r->SetMarkerSize(0.75);
    grW2r->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
    grW2r->GetXaxis()->SetTitleOffset(1.2);
    grW2r->GetXaxis()->CenterTitle();
    grW2r->GetYaxis()->SetTitle("Residuals [keV]");
    grW2r->GetYaxis()->SetTitleOffset(1.2);
    grW2r->GetYaxis()->CenterTitle();
    grW2r->GetXaxis()->SetLimits(0.0,2400.);
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

  else { linCurves << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
    oFileW2 << "PMT NOT USABLE";}
  oFileW2.close();

  // West 3

  sprintf(temp,"../residuals/residuals_West_runPeriod_%i_PMTW3.dat",calibrationPeriod);
  ofstream oFileW3(temp);
  vector <Double_t> fitEQ_W3(runW3.size(),0);

  if (runW3.size()>0 && (std::find(nameW3.begin(),nameW3.end(),"Ce")!=nameW3.end() && std::find(nameW3.begin(),nameW3.end(),"Sn")!=nameW3.end() && std::find(nameW3.begin(),nameW3.end(),"Bi1")!=nameW3.end())) {
    //Calculating the average of the lower Ce ADC values and assigning a value to lowFitThreshold
    Double_t sumADC_Ce = 0., sumADC_Bi = 0., sumEQ_Ce=0.; 
    Int_t entries1 = 0, entries2 = 0;;
    for (UInt_t ii=0; ii<runW3.size(); ii++) {
      if (nameW3[ii]=="Ce") {
	sumADC_Ce+=ADCW3[ii];
	sumEQ_Ce+=EQW3[ii];
	entries1++;
      }
      if (nameW3[ii]=="Bi1") {
	sumADC_Bi+=ADCW3[ii];
	entries2++;
      }
    }

    //Calculating the slope of a line from the origin to the mean of the Ce peaks
    Double_t linSlope = (sumEQ_Ce/sumADC_Ce);
    fitADC->FixParameter(3,linSlope);
    fitADC->FixParameter(4,sumADC_Ce/entries1);
    //lowFitThreshold = sumADC_Ce/entries1;
    highFitThreshold = 1.5*sumADC_Bi/entries2;
    //cout << lowFitThreshold << endl;
    //cout << highFitThreshold << endl;

    
    cW3 = new TCanvas("cW3", "cW3");
    cW3->SetLogy(0);

    TGraphErrors *grW3 = new TGraphErrors(runW3.size(),&ADCW3[0],&EQW3[0],err,err);
    grW3->SetTitle("");
    grW3->SetMarkerColor(1);
    grW3->SetLineColor(1);
    grW3->SetMarkerStyle(20);
    grW3->SetMarkerSize(0.75);
    grW3->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
    grW3->GetXaxis()->SetTitleOffset(1.2);
    grW3->GetXaxis()->CenterTitle();
    grW3->GetYaxis()->SetTitle("E_{Q} [keV]");
    grW3->GetYaxis()->SetTitleOffset(1.6);
    grW3->GetYaxis()->CenterTitle();
    grW3->GetXaxis()->SetLimits(0.0,2400.);
    grW3->SetMinimum(0.0);
    grW3->SetMaximum(1000.0);
    grW3->Draw("AP");

    grW3->Fit("fitADC", "","",0.,highFitThreshold);
  
    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);
    lowEnSlope = fitADC->GetParameter(3);
    shiftOfTanH = fitADC->GetParameter(4);
    spreadOfTanH = fitADC->GetParameter(5);
    slopeToOrigin = (offset + slope*lowFitThreshold + quad*lowFitThreshold*lowFitThreshold + cubic*lowFitThreshold*lowFitThreshold*lowFitThreshold)/lowFitThreshold;
    linCurves << offset << " " << slope << " " << quad << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;
    
    Double_t x1_text = 1200;
    Double_t y1_text = 100;

    TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("West PMT 3");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();

    // Calculate residuals in [keV]
    //Double_t fitEQ_W3[num];
    for (int j=0; j<runW3.size(); j++) {
      fitEQ_W3[j]    = fitADC->Eval(ADCW3[j]);//offset + slope*ADCW3[j] + quad*ADCW3[j]*ADCW3[j] + cubic*ADCW3[j]*ADCW3[j]*ADCW3[j];
      if (nameW3[j]=="Ce") {
	ResW3[j] = fitEQ_W3[j] - EQW3[j];
	oFileW3 << "Ce_West" << " " << runW3[j] << " " << ResW3[j] << endl;
	if (sqrt(ResW3[j]*ResW3[j])>0.025*EQW3[j]) cout<< "Ce_West" << " " << "PMT3 " << runW3[j] << " " << ResW3[j] << endl;  
      }
      else if (nameW3[j]=="Sn") {
	ResW3[j] = fitEQ_W3[j] - EQW3[j];
	oFileW3 << "Sn_West" << " " << runW3[j] << " " << ResW3[j] << endl;
	if (sqrt(ResW3[j]*ResW3[j])>0.025*EQW3[j]) cout<< "Sn_West" << " " << "PMT3 " << runW3[j] << " " << ResW3[j] << endl;  
      }
      else if (nameW3[j]=="Bi1") {
	ResW3[j] = fitEQ_W3[j] - EQW3[j];
	oFileW3 << "Bi1_West" << " " << runW3[j] << " " << ResW3[j] << endl;
	if (sqrt(ResW3[j]*ResW3[j])>0.025*EQW3[j]) cout<< "Bi1_West" << " " << "PMT3 " << runW3[j] << " " << ResW3[j] << endl;  
      }
    }

    // West 3 residuals
    cW3r = new TCanvas("cW3r", "cW3r");
    cW3r->SetLogy(0);

    TGraphErrors *grW3r = new TGraphErrors(runW3.size(),&ADCW3[0],&ResW3[0],err,err);
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
    grW3r->GetXaxis()->SetLimits(0.0,2400.);
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

  else { linCurves << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
    oFileW3 << "PMT NOT USABLE";}
  oFileW3.close();

  // West 4

  sprintf(temp,"../residuals/residuals_West_runPeriod_%i_PMTW4.dat",calibrationPeriod);
  ofstream oFileW4(temp);
  vector <Double_t> fitEQ_W4(runW4.size(),0);
  highFitThreshold=0.;
  if (runW4.size()>0 && (std::find(nameW4.begin(),nameW4.end(),"Ce")!=nameW4.end() && std::find(nameW4.begin(),nameW4.end(),"Sn")!=nameW4.end() && std::find(nameW4.begin(),nameW4.end(),"Bi1")!=nameW4.end())) {
    //Calculating the average of the lower Ce ADC values and assigning a value to lowFitThreshold
    Double_t sumADC_Ce = 0., sumADC_Bi = 0., sumEQ_Ce=0.; 
    Int_t entries1 = 0, entries2 = 0;;
    for (UInt_t ii=0; ii<runW4.size(); ii++) {
      if (nameW4[ii]=="Ce") {
	sumADC_Ce+=ADCW4[ii];
	sumEQ_Ce+=EQW4[ii];
	entries1++;
      }
      if (nameW4[ii]=="Bi1") {
	sumADC_Bi+=ADCW4[ii];
	entries2++;
      }
    }

    //Calculating the slope of a line from the origin to the mean of the Ce peaks
    Double_t linSlope = (sumEQ_Ce/sumADC_Ce);
    fitADC->FixParameter(3,linSlope);
    fitADC->FixParameter(4,sumADC_Ce/entries1);
    //lowFitThreshold = sumADC_Ce/entries1;
    highFitThreshold = 1.5*sumADC_Bi/entries2;
    //cout << lowFitThreshold << endl;
    //cout << highFitThreshold << endl;

    cW4 = new TCanvas("cW4", "cW4");
    cW4->SetLogy(0);

    TGraphErrors *grW4 = new TGraphErrors(runW4.size(),&ADCW4[0],&EQW4[0],err,err);
    grW4->SetTitle("");
    grW4->SetMarkerColor(1);
    grW4->SetLineColor(1);
    grW4->SetMarkerStyle(20);
    grW4->SetMarkerSize(0.75);
    grW4->GetXaxis()->SetTitle("g(t)#upoint#eta(x,y)#upointADC [channels]");
    grW4->GetXaxis()->SetTitleOffset(1.2);
    grW4->GetXaxis()->CenterTitle();
    grW4->GetYaxis()->SetTitle("E_{Q} [keV]");
    grW4->GetYaxis()->SetTitleOffset(1.6);
    grW4->GetYaxis()->CenterTitle();
    grW4->GetXaxis()->SetLimits(0.0,4096.0);
    grW4->SetMinimum(0.0);
    grW4->SetMaximum(1000.0);
    grW4->Draw("AP");

    //fitADC->SetRange(0., 3500.);
    grW4->Fit("fitADC", "","",0.,highFitThreshold);

    offset = fitADC->GetParameter(0);
    slope = fitADC->GetParameter(1);
    quad = fitADC->GetParameter(2);
    lowEnSlope = fitADC->GetParameter(3);
    shiftOfTanH = fitADC->GetParameter(4);
    spreadOfTanH = fitADC->GetParameter(5);
    slopeToOrigin = (offset + slope*lowFitThreshold + quad*lowFitThreshold*lowFitThreshold + cubic*lowFitThreshold*lowFitThreshold*lowFitThreshold)/lowFitThreshold;
    linCurves << offset << " " << slope << " " << quad << " " << lowEnSlope << " " << shiftOfTanH << " " << spreadOfTanH << endl;

    linCurves.close();
  
    Double_t x1_text = 1200;
    Double_t y1_text = 100;

    TPaveText *pt1 = new TPaveText(x1_text,y1_text,x1_text,y1_text,"");
    pt1->SetTextSize(0.042);
    pt1->SetTextColor(1);
    pt1->SetTextAlign(12);
    pt1->AddText("West PMT 4");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    pt1->SetFillColor(0);
    pt1->Draw();

    // Calculate residuals in [keV]
    //Double_t fitEQ_W4[num];
    for (int j=0; j<runW4.size(); j++) {
      fitEQ_W4[j]    = fitADC->Eval(ADCW4[j]);//offset + slope*ADCW4[j] + quad*ADCW4[j]*ADCW4[j] + cubic*ADCW4[j]*ADCW4[j]*ADCW4[j];
      if (nameW4[j]=="Ce") {
	ResW4[j] = fitEQ_W4[j] - EQW4[j];
	oFileW4 << "Ce_West" << " " << runW4[j] << " " << ResW4[j] << endl;
	if (sqrt(ResW4[j]*ResW4[j])>0.025*EQW4[j]) cout<< "Ce_West" << " " << "PMT4 " << runW4[j] << " " << ResW4[j] << endl;  
      }
      else if (nameW4[j]=="Sn") {
	ResW4[j] = fitEQ_W4[j] - EQW4[j];
	oFileW4 << "Sn_West" << " " << runW4[j] << " " << ResW4[j] << endl;
	if (sqrt(ResW4[j]*ResW4[j])>0.025*EQW4[j]) cout<< "Sn_West" << " " << "PMT4 " << runW4[j] << " " << ResW4[j] << endl;  
      }
      else if (nameW4[j]=="Bi1"){
	ResW4[j] = fitEQ_W4[j] - EQW4[j];
	oFileW4 << "Bi1_West" << " " << runW4[j] << " " << ResW4[j] << endl;
	if (sqrt(ResW4[j]*ResW4[j])>0.025*EQW4[j]) cout<< "Bi1_West" << " " << "PMT4 " << runW4[j] << " " << ResW4[j] << endl;  
      }
    }

    // West 4 residuals
    cW4r = new TCanvas("cW4r", "cW4r");
    cW4r->SetLogy(0);

    TGraphErrors *grW4r = new TGraphErrors(runW4.size(),&ADCW4[0],&ResW4[0],err,err);
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
    linCurves << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << " " << 0. << endl;
    oFileW4 << "PMT NOT USABLE";
  }
  oFileW4.close();
  linCurves.close();


  ///////////////////////////////////////////////////////////////////////
  // Calculating the weighted mean of the total energy of the West side  //
  ///////////////////////////////////////////////////////////////////////

  sprintf(temp,"../residuals/residuals_West_runPeriod_%i.dat",calibrationPeriod);
  ofstream oFileW(temp);

  cout << "CALCULATING WEST RESIDUALS" << endl;
  vector <Double_t> Etrue_West(num,0);
  vector <Double_t> x_West(num,0);
  vector <Double_t> res_West(num,0);

  std::vector<Int_t>::iterator W1 = runW1.begin();
  std::vector<Double_t>::iterator W1_EQ = EQW1.begin();
  std::vector<Double_t>::iterator W1_fitEQ = fitEQ_W1.begin();
  std::vector<Int_t>::iterator W2 = runW2.begin();
  std::vector<Double_t>::iterator W2_EQ = EQW2.begin();
  std::vector<Double_t>::iterator W2_fitEQ = fitEQ_W2.begin();
  std::vector<Int_t>::iterator W3 = runW3.begin();
  std::vector<Double_t>::iterator W3_EQ = EQW3.begin();
  std::vector<Double_t>::iterator W3_fitEQ = fitEQ_W3.begin();
  std::vector<Int_t>::iterator W4 = runW4.begin();
  std::vector<Double_t>::iterator W4_EQ = EQW4.begin();
  std::vector<Double_t>::iterator W4_fitEQ = fitEQ_W4.begin();
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
    UInt_t pmtVecLocation = find_vec_location_int(pmtRun,run[j]);

    if (pmtQuality[pmtVecLocation][4] && *W1==run[j] && *nmW1==sourceName[j]) {
      Double_t N = adcW1[j]*nPE_per_channel[pmtVecLocation][4];
      Double_t f = sqrt(N)/N;
      ResidualW1 = *RW1;
      Energy1 = W1_fitEQ;
      std::advance(RW1,1);
      //Energy1 = EQ2Etrue[4][0]+EQ2Etrue[4][1]*(*W1_fitEQ)+EQ2Etrue[4][2]*(*W1_fitEQ)*(*W1_fitEQ);
      weight1 = 1/(Energy1*Energy1*f*f);
      std::advance(W1,1);
      std::advance(W1_EQ,1);
      std::advance(W1_fitEQ,1);
      std::advance(nmW1,1);
    }
    else {weight1=0.; Energy1=0.; ResidualW1=0.;}

    if (pmtQuality[pmtVecLocation][5] && *W2==run[j] && *nmW2==sourceName[j]) {
      Double_t N = adcW2[j]*nPE_per_channel[pmtVecLocation][5];
      Double_t f = sqrt(N)/N;
      ResidualW2 = *RW2;
      Energy2 = W2_fitEQ;
      std::advance(RW2,1);
      //Energy2 = EQ2Etrue[5][0]+EQ2Etrue[5][1]*(*W2_fitEQ)+EQ2Etrue[5][2]*(*W2_fitEQ)*(*W2_fitEQ);
      weight2 = 1/(Energy2*Energy2*f*f);
      std::advance(W2,1);
      std::advance(W2_EQ,1);
      std::advance(W2_fitEQ,1);
      std::advance(nmW2,1);
    }
    else {weight2=0.; Energy2=0.; ResidualW2=0.;}

    if (pmtQuality[pmtVecLocation][6] && *W3==run[j] && *nmW3==sourceName[j]) {
      Double_t N = adcW3[j]*nPE_per_channel[pmtVecLocation][6];
      Double_t f = sqrt(N)/N;
      ResidualW3 = *RW3;
      Energy3 = W3_fitEQ;
      std::advance(RW3,1);
      //Energy3 = EQ2Etrue[6][0]+EQ2Etrue[6][1]*(*W3_fitEQ)+EQ2Etrue[6][2]*(*W3_fitEQ)*(*W3_fitEQ);
      weight3 = 1/(Energy3*Energy3*f*f);
      std::advance(W3,1);
      std::advance(W3_EQ,1);
      std::advance(W3_fitEQ,1);
      std::advance(nmW3,1);
    }
    else {weight3=0.; Energy3=0.; ResidualW3=0.;}

    if (pmtQuality[pmtVecLocation][7] && *W4==run[j] && *nmW4==sourceName[j]) {
      Double_t N = adcW4[j]*nPE_per_channel[pmtVecLocation][7];
      Double_t f = sqrt(N)/N;
      ResidualW4 = *RW4;
      Energy4 = W4_fitEQ;
      std::advance(RW4,1);
      //Energy4 = EQ2Etrue[7][0]+EQ2Etrue[7][1]*(*W4_fitEQ)+EQ2Etrue[7][2]*(*W4_fitEQ)*(*W4_fitEQ);
      weight4 = 1/(Energy4*Energy4*f*f);
      std::advance(W4,1);
      std::advance(W4_EQ,1);
      std::advance(W4_fitEQ,1);
      std::advance(nmW4,1);
    }
    else {weight4=0.; Energy4=0.; ResidualW4=0.;}

    Etrue_West[j] = (weight1*Energy1+weight2*Energy2+weight3*Energy3+weight4*Energy4)/(weight1+weight2+weight3+weight4);
    res_West[j] = (weight1*ResidualW1+weight2*ResidualW2+weight3*ResidualW3+weight4*ResidualW4)/(weight1+weight2+weight3+weight4);
    if (sourceName[j]=="Ce") {
      //res_West[j] = Etrue_West[j] - peakCe;
      x_West[j] = peakCe_EQ;
      oFileW << "Ce_West" << " " <<  run[j] << " " << res_West[j] << endl;
    }
    else if (sourceName[j]=="Sn") {
      //res_West[j] = Etrue_West[j] - peakSn;
      x_West[j] = peakSn_EQ;
      oFileW << "Sn_West" << " " <<  run[j] << " " << res_West[j] << endl;
      /*cout << Energy1 << " " << weight1 << " "
	   << Energy2 << " " << weight2 << " "
	   << Energy3 << " " << weight3 << " " 
	   << Energy4 << " " << weight4 << endl;*/
    }
    else if (sourceName[j]=="Bi2") {
      //res_West[j] = Etrue_West[j] - peakBiLow;
      x_West[j] = peakBiLow_EQ;
    }
    else if (sourceName[j]=="Bi1") {
      //res_West[j] = Etrue_West[j] - peakBiHigh;
      x_West[j] = peakBiHigh_EQ;
      oFileW << "Bi1_West" << " " <<  run[j] << " " << res_West[j] << endl;
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
 
