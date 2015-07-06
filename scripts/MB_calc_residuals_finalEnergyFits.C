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
  // Source peaks from simulation
  Double_t peakCe = 98.2;
  Double_t peakSn = 331.2;
  Double_t peakBiLow = 443.0;
  Double_t peakBiHigh = 928.0;

  // Read data file
  char temp[500];
  sprintf(temp, "../residuals/source_runs_EnergyPeaks_RunPeriod_%i.dat",calibrationPeriod);
  ifstream filein(temp);

  Int_t i = 0;
  Double_t run[500], EQ[500];
  Double_t eastE[500], westE[500];
  Double_t resE[500], resW[500];
  Double_t err[500];
  while (!filein.eof()) {
    filein >> run[i] >> EQ[i]
           >> eastE[i] >> westE[i];
    if (filein.fail()) break;
    i++;
  }
  Int_t num = i;
  cout << "Number of data points: " << num << endl;

  
  sprintf(temp,"../residuals/residuals_East_EnergyPeaks_runPeriod_%i.dat",calibrationPeriod);
  ofstream oFileE(temp);

  vector <double> res_East(num,0);
  vector <double> x_East(num,0);
  
  for (int j=0; j<num; j++) {
  
    if (EQ[j]==98.2) {
      res_East[j] = eastE[j] - peakCe;
      x_East[j] = peakCe;
      oFileE << "Ce_East" << " " << (int) run[j] << " " << res_East[j] << endl;
    }
    else if (EQ[j]==331.2) {
      res_East[j] = eastE[j] - peakSn;
      x_East[j] = peakSn;
      oFileE << "Sn_East" << " " << (int) run[j] << " " << res_East[j] << endl;
    }
    else if (EQ[j]==443.0) {
      res_East[j] = eastE[j] - peakBiLow;
      x_East[j] = peakBiLow;
      oFileE << "Bi2_East" << " " << (int) run[j] << " " << res_East[j] << endl;
    }
    else if (EQ[j]==928.0) {
      res_East[j] = eastE[j] - peakBiHigh;
      x_East[j] = peakBiHigh;
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
    if (EQ[j]==98.2) {
      res_West[j] = westE[j] - peakCe;
      x_West[j] = peakCe;
      oFileW << "Ce_West" << " " << (int) run[j] << " " << res_West[j] << endl;
    }
    else if (EQ[j]==331.2) {
      res_West[j] = westE[j] - peakSn;
      x_West[j] = peakSn;
      oFileW << "Sn_West" << " " << (int) run[j] << " " << res_West[j] << endl;
    }
    else if (EQ[j]==443.0) {
      res_West[j] = westE[j] - peakBiLow;
      x_West[j] = peakBiLow;
      oFileW << "Bi2_West" << " " << (int) run[j] << " " << res_West[j] << endl;
    }
    else if (EQ[j]==928.0) {
      res_West[j] = westE[j] - peakBiHigh;
      x_West[j] = peakBiHigh;
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
