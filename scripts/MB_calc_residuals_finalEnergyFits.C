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

vector<Int_t> pmtRun;

vector < vector <Int_t> > GetPMTQuality(Int_t runPeriod)
{
  vector<vector<int> > pmtQuality;
  vector <int> pmthold(8,0);
  Char_t temp[500];
  
  sprintf(temp,"../residuals/PMT_runQuality_SrcPeriod_%i.dat",runPeriod);
  ifstream pmt;
  pmt.open(temp);
  Int_t run_hold, EPMT1,EPMT2,EPMT3,EPMT4,WPMT1,WPMT2,WPMT3,WPMT4;
  Int_t numRuns=0;
  while (pmt >> run_hold >> pmthold[0] >> pmthold[1] >> pmthold[2]
	 >> pmthold[3] >> pmthold[4] >> pmthold[5]
	 >> pmthold[6] >> pmthold[7]) {
    pmtRun.push_back(run_hold);
    pmtQuality.push_back(pmthold);
    
    numRuns++;
    if (pmt.fail()) break;
  }
  pmt.close();
  return pmtQuality;
}

void MB_calc_residuals_finalEnergyFits(Int_t runPeriod)
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
  sprintf(temp, "../residuals/source_runs_EnergyPeaks_RunPeriod_%i.dat",calibrationPeriod);
  ifstream filein(temp);

  Int_t i = 0;
  Double_t run[500], Etrue[500];
  string sourceName[500];
  Double_t eastE[500], westE[500];
  Double_t resE[500], resW[500];
  Double_t err[500];
  while (!filein.eof()) {
    filein >> run[i] >> sourceName[i]
           >> eastE[i] >> westE[i];
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
 
  //Calculate a weighted average of the smeared EQ peaks weighted by the alpha values
  vector < Double_t > alphas = GetAlphaValues(calibrationPeriod);
  vector < vector <Int_t> > pmtQuality = GetPMTQuality(calibrationPeriod);
  vector < vector <Double_t > > weightedPeaks(4, vector <Double_t> (2,0.));
  Double_t weight[8]={0.};
  Double_t numerE=0., denomE=0., numerW=0., denomW=0.;
  
  //Averaged Smeared Ce peaks
  for (Int_t i=0; i<8; i++) {
    if (pmtQuality[0][i]) {
      weight[i]=alphas[i]/EQsmeared[i][0];
    }
    else weight[i]=0.;
    if (i<4) { numerE+=weight[i]*EQsmeared[i][0]; denomE+=weight[i];}
    else {numerW+=weight[i]*EQsmeared[i][0]; denomW+=weight[i];} 
  }
  weightedPeaks[0][0] = numerE/denomE;
  weightedPeaks[0][1] = numerW/denomW;
  
  //Averaged Smeared Sn peaks
  numerE=denomE=numerW=denomW=0.;
  for (Int_t i=0; i<8; i++) {
    if (pmtQuality[0][i]) {
      weight[i]=alphas[i]/EQsmeared[i][1];
    }
    else weight[i]=0.;
    if (i<4) { numerE+=weight[i]*EQsmeared[i][1]; denomE+=weight[i];}
    if (i>3) {numerW+=weight[i]*EQsmeared[i][1]; denomW+=weight[i];} 
  }
  weightedPeaks[1][0] = numerE/denomE;
  weightedPeaks[1][1] = numerW/denomW;

  //Averaged Smeared Bi2 peaks
  numerE=denomE=numerW=denomW=0.;
  for (Int_t i=0; i<8; i++) {
    if (pmtQuality[0][i]) {
      weight[i]=alphas[i]/EQsmeared[i][2];
    }
    else weight[i]=0.;
    if (i<4) { numerE+=weight[i]*EQsmeared[i][2]; denomE+=weight[i];}
    else {numerW+=weight[i]*EQsmeared[i][2]; denomW+=weight[i];} 
  }
  weightedPeaks[2][0] = numerE/denomE;
  weightedPeaks[2][1] = numerW/denomW;

  //Averaged Smeared Bi1 peaks
  numerE=denomE=numerW=denomW=0.;
  for (Int_t i=0; i<8; i++) {
    if (pmtQuality[0][i]) {
      weight[i]=alphas[i]/EQsmeared[i][3];
    }
    else weight[i]=0.;
    if (i<4) { numerE+=weight[i]*EQsmeared[i][3]; denomE+=weight[i];}
    else {numerW+=weight[i]*EQsmeared[i][3]; denomW+=weight[i];} 
  }
  weightedPeaks[3][0] = numerE/denomE;
  weightedPeaks[3][1] = numerW/denomW;
  cout << weightedPeaks[1][1] << endl;
  ///////////////////////////////////////////////////////////////////


  sprintf(temp,"../residuals/residuals_East_EnergyPeaks_runPeriod_%i.dat",calibrationPeriod);
  ofstream oFileE(temp);

  vector <double> res_East(num,0);
  vector <double> x_East(num,0);
  
  for (int j=0; j<num; j++) {
  
    if (sourceName[j]=="Ce") {
      res_East[j] = eastE[j] - weightedPeaks[0][0];
      x_East[j] = weightedPeaks[0][0];
      oFileE << "Ce_East" << " " << (int) run[j] << " " << res_East[j] << endl;
    }
    else if (sourceName[j]=="Sn") {
      res_East[j] = eastE[j] - weightedPeaks[1][0];
      x_East[j] = weightedPeaks[1][0];
      oFileE << "Sn_East" << " " << (int) run[j] << " " << res_East[j] << endl;
    }
    else if (sourceName[j]=="Bi2") {
      res_East[j] = eastE[j] - weightedPeaks[2][0];
      x_East[j] = weightedPeaks[2][0];
      oFileE << "Bi2_East" << " " << (int) run[j] << " " << res_East[j] << endl;
    }
    else if (sourceName[j]=="Bi1") {
      res_East[j] = eastE[j] - weightedPeaks[3][0];
      x_East[j] = weightedPeaks[3][0];
      oFileE << "Bi1_East" << " " << (int) run[j] << " " << res_East[j] << endl;
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
  grEr->GetXaxis()->SetTitle("E_{true} [keV]");
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
  Double_t x[2] = {0, 2400};
  Double_t y[2] = {0.0, 0.0};

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
  
  ///////////////////////////////////////////////////////////////////////
  // Calculating the weighted mean of the total energy of the West side  //
  ///////////////////////////////////////////////////////////////////////

  sprintf(temp,"../residuals/residuals_West_EnergyPeaks_runPeriod_%i.dat",calibrationPeriod);
  ofstream oFileW(temp);

  cout << "CALCULATING WEST RESIDUALS" << endl;
  vector <Double_t> x_West(num,0);
  vector <Double_t> res_West(num,0);

  for (int j=0; j<num; j++) {
    if (sourceName[j]=="Ce") {
      res_West[j] = westE[j] - weightedPeaks[0][1];
      x_West[j] = weightedPeaks[0][1];
      oFileW << "Ce_West" << " " << (int) run[j] << " " << res_West[j] << endl;
    }
    else if (sourceName[j]=="Sn") {
      res_West[j] = westE[j] - weightedPeaks[1][1];
      x_West[j] = weightedPeaks[1][1];
      oFileW << "Sn_West" << " " << (int) run[j] << " " << res_West[j] << endl;
    }
    else if (sourceName[j]=="Bi2") {
      res_West[j] = westE[j] - weightedPeaks[2][1];
      x_West[j] = weightedPeaks[2][1];
      oFileW << "Bi2_West" << " " << (int) run[j] << " " << res_West[j] << endl;
    }
    else if (sourceName[j]=="Bi1") {
      res_West[j] = westE[j] - weightedPeaks[3][1];
      x_West[j] = weightedPeaks[3][1];
      oFileW << "Bi1_West" << " " << (int) run[j] << " " << res_West[j] << endl;
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
  grWr->GetXaxis()->SetTitle("E_{true} [keV]");
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
  Double_t x[2] = {0, 2400};
  Double_t y[2] = {0.0, 0.0};

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
