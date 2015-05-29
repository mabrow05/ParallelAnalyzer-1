#include <vector>
#include <algorithm>
#include <map>

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

void MB_calc_residuals(Int_t runPeriod)
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
  // Source peaks from simulation
  Double_t peakCe = 98.2;
  Double_t peakSn = 331.2;
  Double_t peakBiLow = 443.0;
  Double_t peakBiHigh = 928.0;

  // Read data file
  char temp[500];
  sprintf(temp, "../residuals/source_runs_RunPeriod_%i.dat",calibrationPeriod);
  ifstream filein(temp);

  Int_t i = 0;
  Double_t run[500], EQ[500];
  Double_t adcE1[500], adcE2[500], adcE3[500], adcE4[500];
  Double_t adcW1[500], adcW2[500], adcW3[500], adcW4[500];
  Double_t resE1[500], resE2[500], resE3[500], resE4[500];
  Double_t resW1[500], resW2[500], resW3[500], resW4[500];
  Double_t err[500];
  while (!filein.eof()) {
    filein >> run[i] >> EQ[i]
           >> adcE1[i] >> adcE2[i] >> adcE3[i] >> adcE4[i]
           >> adcW1[i] >> adcW2[i] >> adcW3[i] >> adcW4[i];
    if (filein.fail()) break;
    i++;
  }
  Int_t num = i;
  cout << "Number of data points: " << num << endl;

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
  vector<Double_t> ADCE1, ADCE2, ADCE3, ADCE4;
  vector<Double_t> ADCW1, ADCW2, ADCW3, ADCW4;
  vector<Double_t> ResE1, ResE2, ResE3, ResE4;
  vector<Double_t> ResW1, ResW2, ResW3, ResW4;
  // Fill run, adc, and eQ vectors
  for (Int_t i=0; i<num;i++) {
    UInt_t runPos = find_vec_location_int(pmtRun,(int)run[i]);
    cout << "Found run " << (int) run[i] << " in PMT Quality\n";
    if (pmtQuality[runPos][0]) {
      runE1.push_back((int)run[i]);
      EQE1.push_back(EQ[i]);
      ADCE1.push_back(adcE1[i]);
    }
    if (pmtQuality[runPos][1]) {
      runE2.push_back((int)run[i]);
      EQE2.push_back(EQ[i]);
      ADCE2.push_back(adcE2[i]);
    }
    if (pmtQuality[runPos][2]) {
      runE3.push_back((int)run[i]);
      EQE3.push_back(EQ[i]);
      ADCE3.push_back(adcE3[i]);
    }
    if (pmtQuality[runPos][3]) {
      runE4.push_back((int)run[i]);
      EQE4.push_back(EQ[i]);
      ADCE4.push_back(adcE4[i]);
    }
    if (pmtQuality[runPos][4]) {
      runW1.push_back((int)run[i]);
      EQW1.push_back(EQ[i]);
      ADCW1.push_back(adcW1[i]);
    }
    if (pmtQuality[runPos][5]) {
      runW2.push_back((int)run[i]);
      EQW2.push_back(EQ[i]);
      ADCW2.push_back(adcW2[i]);
    }
    if (pmtQuality[runPos][6]) {
      runW3.push_back((int)run[i]);
      EQW3.push_back(EQ[i]);
      ADCW3.push_back(adcW3[i]);
    }
    if (pmtQuality[runPos][7]) {
      runW4.push_back((int)run[i]);
      EQW4.push_back(EQ[i]);
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
  

  //File to hold the linearity curves for each calibration run period
  sprintf(temp,"../linearity_curves/lin_curves_srcCal_Period_%i.dat",calibrationPeriod);
  ofstream linCurves(temp);

  // Fit function
  TF1 *fitADC = new TF1("fitADC", "([0] + [1]*x + [2]*x*x)", 0.0, 2500.0);
  fitADC->SetParameter(0, 0.0);
  fitADC->SetParameter(1, 1.0);
  fitADC->SetParameter(2, 0.0);

  fitADC->SetParLimits(2, -0.0001, 0.0001);

  fitADC->SetNpx(100000);
  fitADC->SetLineColor(2);

  // East 1

  sprintf(temp,"../residuals/residuals_East_runPeriod_%i_PMTE1.dat",calibrationPeriod);
  ofstream oFileE1(temp);
  vector <Double_t> fitEQ_E1(runE1.size(),0);

  if (runE1.size()>0 && (std::find(EQE1.begin(),EQE1.end(),peakCe)!=EQE1.end() && std::find(EQE1.begin(),EQE1.end(),peakSn)!=EQE1.end() && std::find(EQE1.begin(),EQE1.end(),peakBiHigh)!=EQE1.end())) {
 
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
    grE1->GetXaxis()->SetLimits(0.0,2400.0);
    grE1->SetMinimum(0.0);
    grE1->SetMaximum(1000.0);
    grE1->Draw("AP");

    grE1->Fit("fitADC", "R");
    Double_t offsetE1, slopeE1, quadE1;
    offsetE1 = fitADC->GetParameter(0);
    slopeE1 = fitADC->GetParameter(1);
    quadE1 = fitADC->GetParameter(2);

    linCurves << offsetE1 << " " << slopeE1 << " " << quadE1 << endl;
  
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
    
      fitEQ_E1[j]    = offsetE1 + slopeE1*ADCE1[j] + quadE1*ADCE1[j]*ADCE1[j];
      if (EQE1[j]==98.2) {
	ResE1[j] = fitEQ_E1[j] - peakCe;
	oFileE1 << "Ce_East" << " " << runE1[j] << " " << ResE1[j] << endl;
	if (ResE1[j]>0.05*peakCe) cout << "Ce_East" << " " << runE1[j] << " " << ResE1[j] << endl;
	//ResE1[j] = (fitEQ - peakCe)/peakCe * 100.;
      }
      else if (EQE1[j]==331.2) {
	ResE1[j] = fitEQ_E1[j] - peakSn;
	//ResE1[j] = (fitEQ - peakSn)/peakSn * 100.;
	oFileE1 << "Sn_East" << " " << runE1[j] << " " << ResE1[j] << endl;
	if (ResE1[j]>0.05*peakSn) cout << "Sn_East" << " " << runE1[j] << " " << ResE1[j] << endl;
      }
      else if (EQE1[j]==928.0) {
	ResE1[j] = fitEQ_E1[j] - peakBiHigh;
	//ResE1[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE1 << "Bi_East" << " " << runE1[j] << " " << ResE1[j] << endl;
	if (ResE1[j]>0.05*peakBiHigh) cout << "Bi_East" << " " << runE1[j] << " " << ResE1[j] << endl;
 
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

  else { linCurves << 0. << " " << 0. << " " << 0. << endl;
    oFileE1 << "PMT NOT USABLE";}
  oFileE1.close();

  // East 2

  sprintf(temp,"../residuals/residuals_East_runPeriod_%i_PMTE2.dat",calibrationPeriod);
  ofstream oFileE2(temp);
  vector <Double_t> fitEQ_E2(runE2.size(),0);

  if (runE2.size()>0 && (std::find(EQE2.begin(),EQE2.end(),peakCe)!=EQE2.end() && std::find(EQE2.begin(),EQE2.end(),peakSn)!=EQE2.end() && std::find(EQE2.begin(),EQE2.end(),peakBiHigh)!=EQE2.end())) {

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
    grE2->GetXaxis()->SetLimits(0.0,2400.0);
    grE2->SetMinimum(0.0);
    grE2->SetMaximum(1000.0);
    grE2->Draw("AP");

    grE2->Fit("fitADC", "R");
    Double_t offsetE2, slopeE2, quadE2;
    offsetE2 = fitADC->GetParameter(0);
    slopeE2 = fitADC->GetParameter(1);
    quadE2 = fitADC->GetParameter(2);

    linCurves << offsetE2 << " " << slopeE2 << " " << quadE2 << endl;

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
      fitEQ_E2[j]    = offsetE2 + slopeE2*ADCE2[j] + quadE2*ADCE2[j]*ADCE2[j];
      if (EQE2[j]==98.2) {
	ResE2[j] = fitEQ_E2[j] - peakCe;
	//ResE2[j] = (fitEQ - peakCe)/peakCe * 100.;
	oFileE2 << "Ce_East" << " " << runE2[j] << " " << ResE2[j] << endl;
	if (ResE2[j]>0.05*peakCe) cout<< "Ce_East" << " " << runE2[j] << " " << ResE2[j] << endl;  
      }
      else if (EQE2[j]==331.2) {
	ResE2[j] = fitEQ_E2[j] - peakSn;
	//ResE2[j] = (fitEQ - peakSn)/peakSn * 100.;
	oFileE2 << "Sn_East" << " " << runE2[j] << " " << ResE2[j] << endl;
	if (ResE2[j]>0.05*peakSn) cout<< "Sn_East" << " " << runE2[j] << " " << ResE2[j] << endl;  
      }
      else if (EQE2[j]==928.0) {
	ResE2[j] = fitEQ_E2[j] - peakBiHigh;
	//ResE2[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE2 << "Bi_East" << " " << runE2[j] << " " << ResE2[j] << endl;
	if (ResE2[j]>0.05*peakCe) cout<< "Bi_East" << " " << runE2[j] << " " << ResE2[j] << endl;  
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

  else { linCurves << 0. << " " << 0. << " " << 0. << endl;
    oFileE2 << "PMT NOT USABLE";}
  oFileE2.close();

  // East 3

  sprintf(temp,"../residuals/residuals_East_runPeriod_%i_PMTE3.dat",calibrationPeriod);
  ofstream oFileE3(temp);
  vector <Double_t> fitEQ_E3(runE3.size(),0);

  if (runE3.size()>0 && (std::find(EQE3.begin(),EQE3.end(),peakCe)!=EQE3.end() && std::find(EQE3.begin(),EQE3.end(),peakSn)!=EQE3.end() && std::find(EQE3.begin(),EQE3.end(),peakBiHigh)!=EQE3.end())) {

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

    grE3->Fit("fitADC", "R");
    Double_t offsetE3, slopeE3, quadE3;
    offsetE3 = fitADC->GetParameter(0);
    slopeE3 = fitADC->GetParameter(1);
    quadE3 = fitADC->GetParameter(2);

    linCurves << offsetE3 << " " << slopeE3 << " " << quadE3 << endl;

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
      fitEQ_E3[j]    = offsetE3 + slopeE3*ADCE3[j] + quadE3*ADCE3[j]*ADCE3[j];
      if (EQE3[j]==98.2) {
	ResE3[j] = fitEQ_E3[j] - peakCe;
	//ResE3[j] = (fitEQ - peakCe)/peakCe * 100.;
	oFileE3 << "Ce_East" << " " << runE3[j] << " " << ResE3[j] << endl;
	if (ResE3[j]>0.05*peakCe) cout<< "Ce_East" << " " << runE3[j] << " " << ResE3[j] << endl;  
      }
      else if (EQE3[j]==331.2) {
	ResE3[j] = fitEQ_E3[j] - peakSn;
	//ResE3[j] = (fitEQ - peakSn)/peakSn * 100.;
	oFileE3 << "Sn_East" << " " << runE3[j] << " " << ResE3[j] << endl;
	if (ResE3[j]>0.05*peakSn) cout << "Sn_East" << " " << runE3[j] << " " << ResE3[j] << endl;  
      }
      else if (EQE3[j]==928.0) {
	ResE3[j] = fitEQ_E3[j] - peakBiHigh;
	//ResE3[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE3 << "Bi_East" << " " << runE3[j] << " " << ResE3[j] << endl;
	if (ResE3[j]>0.05*peakBiHigh) cout << "Bi_East" << " " << runE3[j] << " " << ResE3[j] << endl;  
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

  else { linCurves << 0. << " " << 0. << " " << 0. << endl;
    oFileE3 << "PMT NOT USABLE";}
  oFileE3.close();

  // East 4

  sprintf(temp,"../residuals/residuals_East_runPeriod_%i_PMTE4.dat",calibrationPeriod);
  ofstream oFileE4(temp);
  vector <Double_t> fitEQ_E4(runE4.size(),0);

  if (runE4.size()>0 && (std::find(EQE4.begin(),EQE4.end(),peakCe)!=EQE4.end() && std::find(EQE4.begin(),EQE4.end(),peakSn)!=EQE4.end() && std::find(EQE4.begin(),EQE4.end(),peakBiHigh)!=EQE4.end())) {

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
    grE4->GetXaxis()->SetLimits(0.0,2400.0);
    grE4->SetMinimum(0.0);
    grE4->SetMaximum(1000.0);
    grE4->Draw("AP");

    grE4->Fit("fitADC", "R");
    Double_t offsetE4, slopeE4, quadE4;
    offsetE4 = fitADC->GetParameter(0);
    slopeE4 = fitADC->GetParameter(1);
    quadE4 = fitADC->GetParameter(2);

    linCurves << offsetE4 << " " << slopeE4 << " " << quadE4 << endl;

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
      fitEQ_E4[j]    = offsetE4 + slopeE4*ADCE4[j] + quadE4*ADCE4[j]*ADCE4[j];
      if (EQE4[j]==98.2) {
	ResE4[j] = fitEQ_E4[j] - peakCe;
	//ResE4[j] = (fitEQ - peakCe)/peakCe * 100.;
	oFileE4 << "Ce_East" << " " << runE4[j] << " " << ResE4[j] << endl;
	if (ResE4[j]>0.05*peakCe) cout << "Ce_East" << " " << runE4[j] << " " << ResE4[j] << endl;  
      }
      else if (EQE4[j]==331.2) {
	ResE4[j] = fitEQ_E4[j] - peakSn;
	//ResE4[j] = (fitEQ - peakSn)/peakSn * 100.;
	oFileE4 << "Sn_East" << " " << runE4[j] << " " << ResE4[j] << endl;
	if (ResE4[j]>0.05*peakSn) cout<< "Sn_East" << " " << runE4[j] << " " << ResE4[j] << endl;  
      }
      else if (EQE4[j]==928.0) {
	ResE4[j] = fitEQ_E4[j] - peakBiHigh;
	//ResE4[j] = (fitEQ - peakBiHigh)/peakBiHigh * 100.;
	oFileE4 << "Bi_East" << " " << runE4[j] << " " << ResE4[j] << endl;
	if (ResE4[j]>0.05*peakBiHigh) cout<< "Bi_East" << " " << runE4[j] << " " << ResE4[j] << endl;  
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

  else { linCurves << 0. << " " << 0. << " " << 0. << endl;
    oFileE4 << "PMT NOT USABLE";}
  oFileE4.close();

  ///////////////////////////////////////////////////////////////////////
  // Calcuting the weighted mean of the total energy of the East side  //
  ///////////////////////////////////////////////////////////////////////

  sprintf(temp,"../residuals/residuals_East_runPeriod_%i.dat",calibrationPeriod);
  ofstream oFileE(temp);

  cout << "CALCULATING EAST RESIDUALS" << endl;
  vector <Double_t> EQ_East(num,0);
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

  
  for (int j=0;j<num;j++) {
    Double_t weight1, weight2, weight3, weight4; //these will hold the 4 weights
    Double_t Energy1, Energy2, Energy3, Energy4; //These will hold the energies for each pmt
    UInt_t pmtVecLocation = find_vec_location_int(pmtRun,(int)run[j]);

    if (pmtQuality[pmtVecLocation][0] && *E1==(int)run[j] && *E1_EQ==EQ[j]) {
      Double_t N = adcE1[j]*nPE_per_channel[pmtVecLocation][0];
      Double_t f = sqrt(N)/N;
      Energy1 = *E1_fitEQ;
      weight1 = 1/(Energy1*Energy1*f*f);
      std::advance(E1,1);
      std::advance(E1_EQ,1);
      std::advance(E1_fitEQ,1);
    }
    else {weight1=0.; Energy1=0.;}

    if (pmtQuality[pmtVecLocation][1] && *E2==(int)run[j] && *E2_EQ==EQ[j]) {
      Double_t N = adcE2[j]*nPE_per_channel[pmtVecLocation][1];
      Double_t f = sqrt(N)/N;
      Energy2 = *E2_fitEQ;
      weight2 = 1/(Energy2*Energy2*f*f);
      std::advance(E2,1);
      std::advance(E2_EQ,1);
      std::advance(E2_fitEQ,1);
    }
    else {weight2=0.; Energy2=0.;}

    if (pmtQuality[pmtVecLocation][2] && *E3==(int)run[j] && *E3_EQ==EQ[j]) {
      Double_t N = adcE3[j]*nPE_per_channel[pmtVecLocation][2];
      Double_t f = sqrt(N)/N;
      Energy3 = *E3_fitEQ;
      weight3 = 1/(Energy3*Energy3*f*f);
      std::advance(E3,1);
      std::advance(E3_EQ,1);
      std::advance(E3_fitEQ,1);
    }
    else {weight3=0.; Energy3=0.;}

    if (pmtQuality[pmtVecLocation][3] && *E4==(int)run[j] && *E4_EQ==EQ[j]) {
      Double_t N = adcE4[j]*nPE_per_channel[pmtVecLocation][3];
      Double_t f = sqrt(N)/N;
      Energy4 = *E4_fitEQ;
      weight4 = 1/(Energy4*Energy4*f*f);
      std::advance(E4,1);
      std::advance(E4_EQ,1);
      std::advance(E4_fitEQ,1);
    }
    else {weight4=0.; Energy4=0.;}

    EQ_East[j] = (weight1*Energy1+weight2*Energy2+weight3*Energy3+weight4*Energy4)/(weight1+weight2+weight3+weight4);

    if (EQ[j]==98.2) {
      res_East[j] = EQ_East[j] - peakCe;
      x_East[j] = peakCe;
      oFileE << "Ce_East" << " " << (int) run[j] << " " << res_East[j] << endl;
    }
    else if (EQ[j]==331.2) {
      res_East[j] = EQ_East[j] - peakSn;
      x_East[j] = peakSn;
      oFileE << "Sn_East" << " " << (int) run[j] << " " << res_East[j] << endl;
      cout << Energy1 << " " << weight1 << " "
	   << Energy2 << " " << weight2 << " "
	   << Energy3 << " " << weight3 << " " 
	   << Energy4 << " " << weight4 << endl;
    }
    else if (EQ[j]==443.0) {
      res_East[j] = EQ_East[j] - peakBiLow;
      x_East[j] = peakBiLow;
    }
    else if (EQ[j]==928.0) {
      res_East[j] = EQ_East[j] - peakBiHigh;
      x_East[j] = peakBiHigh;
      oFileE << "Bi_East" << " " << (int) run[j] << " " << res_East[j] << endl;
   }

  }
  oFileE.close();
  
  // East Average residuals
  cEr = new TCanvas("cEr", "cEr");
  cEr->SetLogy(0);

  TGraphErrors *grEr = new TGraphErrors(num,&x_East[0],&res_East[0],err,err);
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
  ////////////////////////////////////////////////////////////////////

  // West 1

  sprintf(temp,"../residuals/residuals_West_runPeriod_%i_PMTW1.dat",calibrationPeriod);
  ofstream oFileW1(temp);
  vector <Double_t> fitEQ_W1(runW1.size(),0);

  if (runW1.size()>0 && (std::find(EQW1.begin(),EQW1.end(),peakCe)!=EQW1.end() && std::find(EQW1.begin(),EQW1.end(),peakSn)!=EQW1.end() && std::find(EQW1.begin(),EQW1.end(),peakBiHigh)!=EQW1.end())) {

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
    grW1->GetXaxis()->SetLimits(0.0,2400.0);
    grW1->SetMinimum(0.0);
    grW1->SetMaximum(1000.0);
    grW1->Draw("AP");

    grW1->Fit("fitADC", "R");
    Double_t offsetW1, slopeW1, quadW1;
    offsetW1 = fitADC->GetParameter(0);
    slopeW1 = fitADC->GetParameter(1);
    quadW1 = fitADC->GetParameter(2);

    linCurves << offsetW1 << " " << slopeW1 << " " << quadW1 << endl;

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
      fitEQ_W1[j]    = offsetW1 + slopeW1*ADCW1[j] + quadW1*ADCW1[j]*ADCW1[j];
      if (EQW1[j]==98.2) {
	ResW1[j] = fitEQ_W1[j] - peakCe;
	oFileW1 << "Ce_West" << " " << runW1[j] << " " << ResW1[j] << endl;
	if (ResW1[j]>0.05*peakCe) cout<< "Ce_West" << " " << runW1[j] << " " << ResW1[j] << endl;  
      }
      else if (EQW1[j]==331.2) {
	ResW1[j] = fitEQ_W1[j] - peakSn;
	oFileW1 << "Sn_West" << " " <<  runW1[j] << " " << ResW1[j] << endl;
	if (ResW1[j]>0.05*peakSn) cout<< "Sn_West" << " " << runW1[j] << " " << ResW1[j] << endl;  
      }
      else if (EQW1[j]==928.0) {
	ResW1[j] = fitEQ_W1[j] - peakBiHigh;
	oFileW1 << "Bi_West" << " " << runW1[j] << " " << ResW1[j] << endl;
	if (ResW1[j]>0.05*peakBiHigh) cout<< "Bi_West" << " " << runW1[j] << " " << ResW1[j] << endl;  
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

  else { linCurves << 0. << " " << 0. << " " << 0. << endl;
    oFileW1 << "PMT NOT USABLE";}
  oFileW1.close();

  // West 2

  sprintf(temp,"../residuals/residuals_West_runPeriod_%i_PMTW2.dat",calibrationPeriod);
  ofstream oFileW2(temp);
  vector <Double_t> fitEQ_W2(runW2.size(),0);

  if (runW2.size()>0 && (std::find(EQW2.begin(),EQW2.end(),peakCe)!=EQW2.end() && std::find(EQW2.begin(),EQW2.end(),peakSn)!=EQW2.end() && std::find(EQW2.begin(),EQW2.end(),peakBiHigh)!=EQW2.end())) {

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
    grW2->GetXaxis()->SetLimits(0.0,2400.0);
    grW2->SetMinimum(0.0);
    grW2->SetMaximum(1000.0);
    grW2->Draw("AP");

    grW2->Fit("fitADC", "R");
    Double_t offsetW2, slopeW2, quadW2;
    offsetW2 = fitADC->GetParameter(0);
    slopeW2 = fitADC->GetParameter(1);
    quadW2 = fitADC->GetParameter(2);

    linCurves << offsetW2 << " " << slopeW2 << " " << quadW2 << endl;

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
      fitEQ_W2[j]    = offsetW2 + slopeW2*ADCW2[j] + quadW2*ADCW2[j]*ADCW2[j];
      if (EQW2[j]==98.2) {
	ResW2[j] = fitEQ_W2[j] - peakCe;
	oFileW2 << "Ce_West" << " " << runW2[j] << " " << ResW2[j] << endl;
	if (ResW2[j]>0.05*peakCe) cout<< "Ce_West" << " " << runW2[j] << " " << ResW2[j] << endl;  
      }
      else if (EQW2[j]==331.2) {
	ResW2[j] = fitEQ_W2[j] - peakSn;
	oFileW2 << "Sn_West" << " " << runW2[j] << " " << ResW2[j] << endl;
	if (ResW2[j]>0.05*peakSn) cout<< "Sn_West" << " " << runW2[j] << " " << ResW2[j] << endl;  
      }
      else if (EQW2[j]==928.0) {
	ResW2[j] = fitEQ_W2[j] - peakBiHigh;
	oFileW2 << "Bi_West" << " " << runW2[j] << " " << ResW2[j] << endl;
	if (ResW2[j]>0.05*peakBiHigh) cout<< "Bi_West" << " " << runW2[j] << " " << ResW2[j] << endl;  
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

  else { linCurves << 0. << " " << 0. << " " << 0. << endl;
    oFileW2 << "PMT NOT USABLE";}
  oFileW2.close();

  // West 3

  sprintf(temp,"../residuals/residuals_West_runPeriod_%i_PMTW3.dat",calibrationPeriod);
  ofstream oFileW3(temp);
  vector <Double_t> fitEQ_W3(runW3.size(),0);

  if (runW3.size()>0 && (std::find(EQW3.begin(),EQW3.end(),peakCe)!=EQW3.end() && std::find(EQW3.begin(),EQW3.end(),peakSn)!=EQW3.end() && std::find(EQW3.begin(),EQW3.end(),peakBiHigh)!=EQW3.end())) {

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
    grW3->GetXaxis()->SetLimits(0.0,2400.0);
    grW3->SetMinimum(0.0);
    grW3->SetMaximum(1000.0);
    grW3->Draw("AP");

    grW3->Fit("fitADC", "R");
    Double_t offsetW3, slopeW3, quadW3;
    offsetW3 = fitADC->GetParameter(0);
    slopeW3 = fitADC->GetParameter(1);
    quadW3 = fitADC->GetParameter(2);

    linCurves << offsetW3 << " " << slopeW3 << " " << quadW3 << endl;
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
      fitEQ_W3[j]    = offsetW3 + slopeW3*ADCW3[j] + quadW3*ADCW3[j]*ADCW3[j];
      if (EQW3[j]==98.2) {
	ResW3[j] = fitEQ_W3[j] - peakCe;
	oFileW3 << "Ce_West" << " " << runW3[j] << " " << ResW3[j] << endl;
	if (ResW3[j]>0.05*peakCe) cout<< "Ce_West" << " " << runW3[j] << " " << ResW3[j] << endl;  
      }
      else if (EQW3[j]==331.2) {
	ResW3[j] = fitEQ_W3[j] - peakSn;
	oFileW3 << "Sn_West" << " " << runW3[j] << " " << ResW3[j] << endl;
	if (ResW3[j]>0.05*peakSn) cout<< "Sn_West" << " " << runW3[j] << " " << ResW3[j] << endl;  
      }
      else if (EQW3[j]==928.0) {
	ResW3[j] = fitEQ_W3[j] - peakBiHigh;
	oFileW3 << "Bi_West" << " " << runW3[j] << " " << ResW3[j] << endl;
	if (ResW3[j]>0.05*peakBiHigh) cout<< "Bi_West" << " " << runW3[j] << " " << ResW3[j] << endl;  
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

  else { linCurves << 0. << " " << 0. << " " << 0. << endl;
    oFileW3 << "PMT NOT USABLE";}
  oFileW3.close();

  // West 4

  sprintf(temp,"../residuals/residuals_West_runPeriod_%i_PMTW4.dat",calibrationPeriod);
  ofstream oFileW4(temp);
  vector <Double_t> fitEQ_W4(runW4.size(),0);

  if (runW4.size()>0 && (std::find(EQW4.begin(),EQW4.end(),peakCe)!=EQW4.end() && std::find(EQW4.begin(),EQW4.end(),peakSn)!=EQW4.end() && std::find(EQW4.begin(),EQW4.end(),peakBiHigh)!=EQW4.end())) {

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

    fitADC->SetRange(0., 3500.);
    grW4->Fit("fitADC", "R");
    Double_t offsetW4, slopeW4, quadW4;
    offsetW4 = fitADC->GetParameter(0);
    slopeW4 = fitADC->GetParameter(1);
    quadW4 = fitADC->GetParameter(2);

    linCurves << offsetW4 << " " << slopeW4 << " " << quadW4;
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
      fitEQ_W4[j]    = offsetW4 + slopeW4*ADCW4[j] + quadW4[j]*ADCW4[j]*ADCW4[j];
      if (EQW4[j]==98.2) {
	ResW4[j] = fitEQ_W4[j] - peakCe;
	oFileW4 << "Ce_West" << " " << runW4[j] << " " << ResW4[j] << endl;
	if (ResW4[j]>0.05*peakCe) cout<< "Ce_West" << " " << runW4[j] << " " << ResW4[j] << endl;  
      }
      else if (EQW4[j]==331.2) {
	ResW4[j] = fitEQ_W4[j] - peakSn;
	oFileW4 << "Sn_West" << " " << runW4[j] << " " << ResW4[j] << endl;
	if (ResW4[j]>0.05*peakSn) cout<< "Sn_West" << " " << runW4[j] << " " << ResW4[j] << endl;  
      }
      else if (EQW4[j]==928.0){
	ResW4[j] = fitEQ_W4[j] - peakBiHigh;
	oFileW4 << "Bi_West" << " " << runW4[j] << " " << ResW4[j] << endl;
	if (ResW4[j]>0.05*peakBiHigh) cout<< "Bi_West" << " " << runW4[j] << " " << ResW4[j] << endl;  
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
    linCurves << 0. << " " << 0. << " " << 0. << endl;
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
  vector <Double_t> EQ_West(num,0);
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

  
  for (int j=0;j<num;j++) {
    Double_t weight1, weight2, weight3, weight4; //these will hold the 4 weights
    Double_t Energy1, Energy2, Energy3, Energy4; //These will hold the energies for each pmt
    UInt_t pmtVecLocation = find_vec_location_int(pmtRun,(int)run[j]);

    if (pmtQuality[pmtVecLocation][4] && *W1==(int)run[j] && *W1_EQ==EQ[j]) {
      Double_t N = adcW1[j]*nPE_per_channel[pmtVecLocation][4];
      Double_t f = sqrt(N)/N;
      Energy1 = *W1_fitEQ;
      weight1 = 1/(Energy1*Energy1*f*f);
      std::advance(W1,1);
      std::advance(W1_EQ,1);
      std::advance(W1_fitEQ,1);
    }
    else {weight1=0.; Energy1=0.;}

    if (pmtQuality[pmtVecLocation][5] && *W2==(int)run[j] && *W2_EQ==EQ[j]) {
      Double_t N = adcW2[j]*nPE_per_channel[pmtVecLocation][5];
      Double_t f = sqrt(N)/N;
      Energy2 = *W2_fitEQ;
      weight2 = 1/(Energy2*Energy2*f*f);
      std::advance(W2,1);
      std::advance(W2_EQ,1);
      std::advance(W2_fitEQ,1);
    }
    else {weight2=0.; Energy2=0.;}

    if (pmtQuality[pmtVecLocation][6] && *W3==(int)run[j] && *W3_EQ==EQ[j]) {
      Double_t N = adcW3[j]*nPE_per_channel[pmtVecLocation][6];
      Double_t f = sqrt(N)/N;
      Energy3 = *W3_fitEQ;
      weight3 = 1/(Energy3*Energy3*f*f);
      std::advance(W3,1);
      std::advance(W3_EQ,1);
      std::advance(W3_fitEQ,1);
    }
    else {weight3=0.; Energy3=0.;}

    if (pmtQuality[pmtVecLocation][7] && *W4==(int)run[j] && *W4_EQ==EQ[j]) {
      Double_t N = adcW4[j]*nPE_per_channel[pmtVecLocation][7];
      Double_t f = sqrt(N)/N;
      Energy4 = *W4_fitEQ;
      weight4 = 1/(Energy4*Energy4*f*f);
      std::advance(W4,1);
      std::advance(W4_EQ,1);
      std::advance(W4_fitEQ,1);
    }
    else {weight4=0.; Energy4=0.;}

    EQ_West[j] = (weight1*Energy1+weight2*Energy2+weight3*Energy3+weight4*Energy4)/(weight1+weight2+weight3+weight4);

    if (EQ[j]==98.2) {
      res_West[j] = EQ_West[j] - peakCe;
      x_West[j] = peakCe;
      oFileW << "Ce_West" << " " << (int) run[j] << " " << res_West[j] << endl;
    }
    else if (EQ[j]==331.2) {
      res_West[j] = EQ_West[j] - peakSn;
      x_West[j] = peakSn;
      oFileW << "Sn_West" << " " << (int) run[j] << " " << res_West[j] << endl;
      cout << Energy1 << " " << weight1 << " "
	   << Energy2 << " " << weight2 << " "
	   << Energy3 << " " << weight3 << " " 
	   << Energy4 << " " << weight4 << endl;
    }
    else if (EQ[j]==443.0) {
      res_West[j] = EQ_West[j] - peakBiLow;
      x_West[j] = peakBiLow;
    }
    else if (EQ[j]==928.0) {
      res_West[j] = EQ_West[j] - peakBiHigh;
      x_West[j] = peakBiHigh;
      oFileW << "Bi_West" << " " << (int) run[j] << " " << res_West[j] << endl;
   }

  }
  oFileW.close();

  // West Average residuals
  cWr = new TCanvas("cWr", "cWr");
  cWr->SetLogy(0);

  TGraphErrors *grWr = new TGraphErrors(num,&x_West[0],&res_West[0],err,err);
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
  
  //outfile->Write();
  //outfile->Close();

}
