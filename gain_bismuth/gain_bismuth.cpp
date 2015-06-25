#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include <sstream>

// ROOT libraries
#include "TRandom3.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TH1D.h>

#include "fullTreeVariables.h"
#include "MWPCGeometry.h"
#include "pedestals.h"
#include "cuts.h"
#include "basic_reconstruction.h"

using namespace std;

int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Open output ntuple
  char tempOut[500];
  sprintf(tempOut, "%s/gain_bismuth_%s.root",getenv("GAIN_BISMUTH"),argv[1]);
  TFile *fileOut = new TFile(tempOut,"RECREATE");

  // Run number integer
  istringstream ss(argv[1]);
  int runNumber;
  ss >> runNumber;
  cout << "runNumber = " << runNumber << endl;

  // Output histograms
  int nBin = 200;
  TH1F *his[8];
  his[0] = new TH1F("hisE0", "", nBin,0.0,4000.0);
  his[1] = new TH1F("hisE1", "", nBin,0.0,4000.0);
  his[2] = new TH1F("hisE2", "", nBin,0.0,4000.0);
  his[3] = new TH1F("hisE3", "", nBin,0.0,4000.0);
  his[4] = new TH1F("hisW0", "", nBin,0.0,4000.0);
  his[5] = new TH1F("hisW1", "", nBin,0.0,4000.0);
  his[6] = new TH1F("hisW2", "", nBin,0.0,4000.0);
  his[7] = new TH1F("hisW3", "", nBin,0.0,4000.0);

  // Open input ntuple
  char tempIn[500];
  sprintf(tempIn, "%s/replay_pass1_%s.root",getenv("REPLAY_PASS1"), argv[1]);
  TFile *fileIn = new TFile(tempIn, "READ");
  TTree *Tin = (TTree*)(fileIn->Get("pass1"));
  int nEvents = Tin->GetEntries();
  cout << "Processing " << argv[1] << " ... " << endl;
  cout << "... nEvents = " << nEvents << endl;

  // Loop over events
  for (int i=0; i<nEvents; i++) {
    Tin->GetEvent(i);

    // Input variables
    Tin->SetBranchAddress("pmt0", &pmt[0]);
    Tin->SetBranchAddress("pmt1", &pmt[1]);
    Tin->SetBranchAddress("pmt2", &pmt[2]);
    Tin->SetBranchAddress("pmt3", &pmt[3]);
    Tin->SetBranchAddress("pmt4", &pmt[4]);
    Tin->SetBranchAddress("pmt5", &pmt[5]);
    Tin->SetBranchAddress("pmt6", &pmt[6]);
    Tin->SetBranchAddress("pmt7", &pmt[7]);

    Tin->SetBranchAddress("PID", &PID);
    Tin->SetBranchAddress("type", &type);
    Tin->SetBranchAddress("side", &side);

    // Select Bi pulser events
    if (PID != 4) continue;

    // Fill PMT histograms
    his[0]->Fill(pmt[0]);
    his[1]->Fill(pmt[1]);
    his[2]->Fill(pmt[2]);
    his[3]->Fill(pmt[3]);
    his[4]->Fill(pmt[4]);
    his[5]->Fill(pmt[5]);
    his[6]->Fill(pmt[6]);
    his[7]->Fill(pmt[7]);

  }

  // Find maximum bin
  double maxBin[8];
  double maxCounts[8];
  for (int i=0; i<8; i++) {
    maxCounts[i] = -1.0;
  }

  double binCenter[8];
  double binCenterMax[8];
  double binCounts[8];
  for (int i=0; i<nBin; i++) {
    for (int j=0; j<8; j++) {
      binCenter[j] = his[j]->GetBinCenter(i+1);
      binCounts[j] = his[j]->GetBinContent(i+1);
      if (binCenter[j] > 500. && binCenter[j] < 3500. && binCounts[j] >= maxCounts[j]) {
        maxCounts[j] = binCounts[j];
        maxBin[j] = i;
        binCenterMax[j] = binCenter[j];
      }
    }
  }

  // Define histogram fit ranges
  double xLow[8], xHigh[8];
  for (int n=0; n<8; n++) {
    for (int i=maxBin[n]; i<nBin; i++) {
      if (his[n]->GetBinContent(i+1) < 0.4*maxCounts[n]) {
        xHigh[n] = his[n]->GetBinCenter(i+1);
        break;
      }
    }
    for (int i=maxBin[n]; i>0; i--) {
      if (his[n]->GetBinContent(i+1) < 0.4*maxCounts[n]) {
        xLow[n] = his[n]->GetBinCenter(i+1);
        break;
      }
    }
  }

  // Fit parameters
  double fitMean[8];
  for (int m=0; m<8; m++) {
    fitMean[m] = -1.0;
  }

  // Fit for East PMT #0
  TF1 *gaussian_fit_E0 = new TF1("gaussian_fit_E0",
                                 "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
                                 xLow[0], xHigh[0]);
  gaussian_fit_E0->SetParameter(0, maxCounts[0]);
  gaussian_fit_E0->SetParameter(1, binCenterMax[0]);
  gaussian_fit_E0->SetParameter(2, 200.);
  gaussian_fit_E0->SetLineColor(2);
  his[0]->Fit("gaussian_fit_E0", "LR");
  fitMean[0] = gaussian_fit_E0->GetParameter(1);
  cout << "fitMean[0] = " << fitMean[0] << endl;

  // Fit for East PMT #1
  TF1 *gaussian_fit_E1 = new TF1("gaussian_fit_E1",
                                 "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
                                 xLow[1], xHigh[1]);
  gaussian_fit_E1->SetParameter(0, maxCounts[1]);
  gaussian_fit_E1->SetParameter(1, binCenterMax[1]);
  gaussian_fit_E1->SetParameter(2, 200.);
  gaussian_fit_E1->SetLineColor(2);
  his[1]->Fit("gaussian_fit_E1", "LR");
  fitMean[1] = gaussian_fit_E1->GetParameter(1);
  cout << "fitMean[1] = " << fitMean[1] << endl;

  // Fit for East PMT #2
  TF1 *gaussian_fit_E2 = new TF1("gaussian_fit_E2",
                                 "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
                                 xLow[2], xHigh[2]);
  gaussian_fit_E2->SetParameter(0, maxCounts[2]);
  gaussian_fit_E2->SetParameter(1, binCenterMax[2]);
  gaussian_fit_E2->SetParameter(2, 200.);
  gaussian_fit_E2->SetLineColor(2);
  his[2]->Fit("gaussian_fit_E2", "LR");
  fitMean[2] = gaussian_fit_E2->GetParameter(1);
  cout << "fitMean[2] = " << fitMean[2] << endl;

  // Fit for East PMT #3
  TF1 *gaussian_fit_E3 = new TF1("gaussian_fit_E3",
                                 "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
                                 xLow[3], xHigh[3]);
  gaussian_fit_E3->SetParameter(0, maxCounts[3]);
  gaussian_fit_E3->SetParameter(1, binCenterMax[3]);
  gaussian_fit_E3->SetParameter(2, 200.);
  gaussian_fit_E3->SetLineColor(2);
  his[3]->Fit("gaussian_fit_E3", "LR");
  fitMean[3] = gaussian_fit_E3->GetParameter(1);
  cout << "fitMean[3] = " << fitMean[3] << endl;

  // Fit for West PMT #0
  TF1 *gaussian_fit_W0 = new TF1("gaussian_fit_W0",
                                 "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
                                 xLow[4], xHigh[4]);
  gaussian_fit_W0->SetParameter(0, maxCounts[4]);
  gaussian_fit_W0->SetParameter(1, binCenterMax[4]);
  gaussian_fit_W0->SetParameter(2, 200.);
  gaussian_fit_W0->SetLineColor(2);
  his[4]->Fit("gaussian_fit_W0", "LR");
  fitMean[4] = gaussian_fit_W0->GetParameter(1);
  cout << "fitMean[4] = " << fitMean[4] << endl;

  // Fit for West PMT #1
  TF1 *gaussian_fit_W1 = new TF1("gaussian_fit_W1",
                                 "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
                                 xLow[5], xHigh[5]);
  gaussian_fit_W1->SetParameter(0, maxCounts[5]);
  gaussian_fit_W1->SetParameter(1, binCenterMax[5]);
  gaussian_fit_W1->SetParameter(2, 200.);
  gaussian_fit_W1->SetLineColor(2);
  his[5]->Fit("gaussian_fit_W1", "LR");
  fitMean[5] = gaussian_fit_W1->GetParameter(1);
  cout << "fitMean[5] = " << fitMean[5] << endl;

  // Fit for West PMT #2
  TF1 *gaussian_fit_W2 = new TF1("gaussian_fit_W2",
                                 "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
                                 xLow[6], xHigh[6]);
  gaussian_fit_W2->SetParameter(0, maxCounts[6]);
  gaussian_fit_W2->SetParameter(1, binCenterMax[6]);
  gaussian_fit_W2->SetParameter(2, 200.);
  gaussian_fit_W2->SetLineColor(2);
  his[6]->Fit("gaussian_fit_W2", "LR");
  fitMean[6] = gaussian_fit_W2->GetParameter(1);
  cout << "fitMean[6] = " << fitMean[6] << endl;

  // Fit for West PMT #3
  TF1 *gaussian_fit_W3 = new TF1("gaussian_fit_W3",
                                 "[0]*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))",
                                 xLow[7], xHigh[7]);
  gaussian_fit_W3->SetParameter(0, maxCounts[7]);
  gaussian_fit_W3->SetParameter(1, binCenterMax[7]);
  gaussian_fit_W3->SetParameter(2, 200.);
  gaussian_fit_W3->SetLineColor(2);
  his[7]->Fit("gaussian_fit_W3", "LR");
  fitMean[7] = gaussian_fit_W3->GetParameter(1);
  cout << "fitMean[7] = " << fitMean[7] << endl;

  // Calculate gain correction factors
  double referenceMean[8];

  // Source Calibration Run Period 1
  if (runNumber >= 16983 && runNumber <= 17297) {
    referenceMean[0] = 2768.38;
    referenceMean[1] = 2911.87;
    referenceMean[2] = 2820.94;
    referenceMean[3] = 3020.11;
    referenceMean[4] = 2870.89;
    referenceMean[5] = -1.;
    referenceMean[6] = 2756.86;
    referenceMean[7] = 2887.47;
  }

  // Source Calibration Run Period 2
  if (runNumber >= 17359 && runNumber <= 17439) {
    referenceMean[0] = 2783.11;
    referenceMean[1] = 2825.62;
    referenceMean[2] = 2665.44;
    referenceMean[3] = 2613.36;
    referenceMean[4] = 715.844;
    referenceMean[5] = 2471.9;
    referenceMean[6] = 2497.24;
    referenceMean[7] = 2974.43;
  }

  // Source Calibration Run Period 3
  if (runNumber >= 17440 && runNumber <= 17734) {
    referenceMean[0] = 2925.49;
    referenceMean[1] = 2860.85;
    referenceMean[2] = 2646.09;
    referenceMean[3] = 2609.76;
    referenceMean[4] = 728.335;
    referenceMean[5] = 2740.75;
    referenceMean[6] = 2591.74;
    referenceMean[7] = 2458.67;
  }

  // Source Calibration Run Period 4
  if (runNumber >= 17735 && runNumber <= 18055) {
    referenceMean[0] = 2913.21;
    referenceMean[1] = 2726.86;
    referenceMean[2] = 2614.41;
    referenceMean[3] = 2574.88;
    referenceMean[4] = 685.297;
    referenceMean[5] = 2710.91;
    referenceMean[6] = 2557.8;
    referenceMean[7] = 2457.16;
  }

  // Source Calibration Run Period 5
  if (runNumber >= 18081 && runNumber <= 18386) {
    referenceMean[0] = 3139.68;
    referenceMean[1] = 2478.71;
    referenceMean[2] = 2635.56;
    referenceMean[3] = 2603.68;
    referenceMean[4] = 2602.94;
    referenceMean[5] = 2118.06;
    referenceMean[6] = 2549.17;
    referenceMean[7] = 1710.95;
  }

  // Source Calibration Run Period 6
  if (runNumber >= 18390 && runNumber <= 18683) {
    referenceMean[0] = 2798.3;
    referenceMean[1] = 2648.22;
    referenceMean[2] = 2643.97;
    referenceMean[3] = 2580.29;
    referenceMean[4] = 2557.37;
    referenceMean[5] = 2850.5;
    referenceMean[6] = 2484.63;
    referenceMean[7] = 1266.61;
  }

  // Source Calibration Run Period 7
  if (runNumber >= 18712 && runNumber <= 18994) {
    referenceMean[0] = 2473.59;
    referenceMean[1] = 2490.29;
    referenceMean[2] = 2635.3;
    referenceMean[3] = 2575.75;
    referenceMean[4] = 2562.27;
    referenceMean[5] = 2609.88;
    referenceMean[6] = 2456.6;
    referenceMean[7] = 1313.06;
  }

  // Source Calibration Run Period 8
  if (runNumber >= 19203 && runNumber <= 19239) {
    referenceMean[0] = 2543.39;
    referenceMean[1] = 2581.96;
    referenceMean[2] = 2697.13;
    referenceMean[3] = 2607.16;
    referenceMean[4] = 2583.02;
    referenceMean[5] = 2649.35;
    referenceMean[6] = 2467.87;
    referenceMean[7] = 2617.73;
  }

  // Source Calibration Run Period 9
  if (runNumber >= 19347 && runNumber <= 19544) {
    referenceMean[0] = 2532.88;
    referenceMean[1] = 2442.47;
    referenceMean[2] = 2685.81;
    referenceMean[3] = 2668.67;
    referenceMean[4] = 2592.06;
    referenceMean[5] = 2628.69;
    referenceMean[6] = 2520.98;
    referenceMean[7] = 781.39; //This PMT is not used for the calibration
  }


  // Source Calibration Run Period 10
  // This calibration isn't used as of now. These reference means are 
  // not correct! If this is to ever be used, we need to look up the 
  // values as fit for the reference run in this period.
  // FOR NOW, USE NEXT RUN PERIOD GAIN (no Ce/Sn/Bi in Run Period 10)
  /*if (runNumber >= 19505 && runNumber <= 19544) {
    referenceMean[0] = 2620.29;
    referenceMean[1] = 2530.83;
    referenceMean[2] = 2689.98;
    referenceMean[3] = 2657.57;
    referenceMean[4] = 2556.1;
    referenceMean[5] = 2627.63;
    referenceMean[6] = 2531.34;
    referenceMean[7] = 797.798;
    }*/

  // Source Calibration Run Period 11
  if (runNumber >= 19583 && runNumber <= 20000) {
    referenceMean[0] = 2620.29;
    referenceMean[1] = 2530.83;
    referenceMean[2] = 2689.98;
    referenceMean[3] = 2657.57;
    referenceMean[4] = 2556.1;
    referenceMean[5] = 2627.63;
    referenceMean[6] = 2531.34;
    referenceMean[7] = 797.798;
  }

  // Source Calibration Run Period 12
  /*if (runNumber >= 19899 && runNumber <= 19966) {
    referenceMean[0] = 2627.55;
    referenceMean[1] = 2540.88;
    referenceMean[2] = 2688.1;
    referenceMean[3] = 2646.31;
    referenceMean[4] = 2560.41;
    referenceMean[5] = 2642.34;
    referenceMean[6] = 2537.37;
    referenceMean[7] = 776.104;
    }*/

  // Calculate gain correction factors
  double gainCorrection[8];
  for (int n=0; n<8; n++) {
    gainCorrection[n] = referenceMean[n] / fitMean[n];
  }

  // Close input ntuple
  fileIn->Close();

  // Write output ntuple
  fileOut->Write();
  fileOut->Close();

  // Write fit results to file
  char tempFit[500];
  sprintf(tempFit, "%s/gain_bismuth_%s.dat", getenv("GAIN_BISMUTH"),argv[1]);
  ofstream outFit(tempFit);

  for (int n=0; n<8; n++) {
    //outFit << gainCorrection[n] << endl;
    outFit << fitMean[n] << " " << gainCorrection[n] << endl;
  }

  return 0;
}
