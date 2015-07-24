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

#include "replay_pass2.h"
#include "replay_pass3.h"
#include "replay_pass4.h"
#include "sourcePeaks.h"

using namespace std;

int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);  

  Int_t nPMTs = 8;
  Int_t nParams = 3; //quadratic fit
  
  // Run number integer
  cout << "Run " << argv[1] << " ..." << endl;
  cout << "... Applying Energy Calibration ..." << endl;
  //istringstream ss(argv[1]);
  int runNumber = atoi(argv[1]);
  //ss >> runNumber;
  int calibrationPeriod = 1;
  // Determine linearity curve to use
  char tempFileLinearityCurve[500];
  if (runNumber <= 17297) {
    calibrationPeriod=1;
    sprintf(tempFileLinearityCurve, "../linearity_curves/lin_curves_srcCal_Period_1.dat");
  }
  else if (runNumber <= 17439) {
    calibrationPeriod=2;
    sprintf(tempFileLinearityCurve, "../linearity_curves/lin_curves_srcCal_Period_2.dat");
  }
  else if (runNumber <= 17734) {
    calibrationPeriod=3;
    sprintf(tempFileLinearityCurve, "../linearity_curves/lin_curves_srcCal_Period_3.dat");
  }
  else if (runNumber <= 17955) {
    calibrationPeriod=4;
    sprintf(tempFileLinearityCurve, "../linearity_curves/lin_curves_srcCal_Period_4.dat");
  }
  else if (runNumber <= 18386) {
    calibrationPeriod=5;
    sprintf(tempFileLinearityCurve, "../linearity_curves/lin_curves_srcCal_Period_5.dat");
  }
  else if (runNumber <= 18683) {
    calibrationPeriod=6;
    sprintf(tempFileLinearityCurve, "../linearity_curves/lin_curves_srcCal_Period_6.dat");
  }
  else if (runNumber <= 18994) {
    calibrationPeriod=7;
    sprintf(tempFileLinearityCurve, "../linearity_curves/lin_curves_srcCal_Period_7.dat");
  }
  else if (runNumber <= 19239) {
    calibrationPeriod=8;
    sprintf(tempFileLinearityCurve, "../linearity_curves/lin_curves_srcCal_Period_8.dat");
  }
  else if (runNumber <= 19544) {
    calibrationPeriod=9;
    sprintf(tempFileLinearityCurve, "../linearity_curves/lin_curves_srcCal_Period_9.dat");
  }
  else if (runNumber <= 20000) {
    calibrationPeriod=11;
    sprintf(tempFileLinearityCurve, "../linearity_curves/lin_curves_srcCal_Period_11.dat");
  }
  cout << "... Reading: " << tempFileLinearityCurve << endl;

   //Setup array to hold linearity curves
  Double_t linearityCurve[nPMTs][nParams];

  // Read linearity curve
  cout << "Reading in linearity curve:\n";
  cout << "p0\tp1\tp2\n";
  ifstream fileLinearityCurve(tempFileLinearityCurve);
  Double_t p0,p1,p2;
  Int_t i=0;
  while (fileLinearityCurve >> p0 >> p1 >> p2) {
    linearityCurve[i][0] = p0;
    linearityCurve[i][1] = p1;
    linearityCurve[i][2] = p2;
    i++;
    cout << p0 << " " << p1 << " " << p2 << endl;
    if (fileLinearityCurve.fail()) break;                       
  }

  //Load the simulated relationship between EQ and Etrue
  vector < vector <double> > EQ2Etrue = EQ2EtrueFit(calibrationPeriod);
  cout << "EQ to Etrue conversion quadratic:\n" << "const:\tlin:\tquad:\n";
  for (unsigned int m=0; m<EQ2Etrue.size(); m++) {
    for (unsigned int mm=0; mm<EQ2Etrue[m].size(); mm++) {
      cout << EQ2Etrue[m][mm] << " ";
    }
    cout << endl;
  }
  
  //Read in PMT quality file
  cout << "Reading in PMT Quality file ...\n";
  Int_t pmtQuality[8];
  Char_t temp[200];
  sprintf(temp,"../residuals/PMT_runQuality_master.dat"); 
  ifstream pmt;
  pmt.open(temp);
  Int_t run_hold;
  while (pmt >> run_hold >> pmtQuality[0] >> pmtQuality[1] >> pmtQuality[2]
	 >> pmtQuality[3] >> pmtQuality[4] >> pmtQuality[5]
	 >> pmtQuality[6] >> pmtQuality[7]) {
    if (run_hold==runNumber) break;
    if (pmt.fail()) break;
  }
  pmt.close();
  if (run_hold!=runNumber) {
    cout << "Run not found in PMT quality file!" << endl;
    exit(0);
  }
  //Also must check that the linearity curves being used had all functioning PMTs... 
  for (int j=0; j<8; j++) {
    if (linearityCurve[j][1]==0.) pmtQuality[j]=0;
  }

  cout << pmtQuality[0] << " " << pmtQuality[1] << " " << pmtQuality[2] 
       << " " << pmtQuality[3] << " " << pmtQuality[4] << " " << pmtQuality[5] 
       << " " << pmtQuality[6] << " " << pmtQuality[7] << endl; 

  //Read in PMT weights
  cout << "Reading in PMT weighting factors...\n";
  // Filling vectors with the weights for each PMT to be used when calculating the weighted average energy
 
  Double_t nPE_per_channel[nPMTs];
  Double_t nPE_hold, sigma_hold, mean_hold, nPE_per_chan_hold;
  ifstream weightFile;
  sprintf(temp,"%s/nPE_weights_%i.dat",getenv("NPE_WEIGHTS"),runNumber);
  weightFile.open(temp);
  int ii=0;
  while (weightFile >> mean_hold >> sigma_hold >> nPE_hold >> nPE_per_chan_hold) {
    nPE_per_channel[ii]=nPE_per_chan_hold;
    ii++;
    cout << ii << " " << nPE_per_channel[ii-1] << endl;
    if (weightFile.fail()) break;
  }
  weightFile.close();

  if (ii!=8) {
    cout << "Didn't read in PMT Weights for all 8 PMTs\n";
    exit(0);
  }

  // Open output ntuple
  char tempOut[500];
  sprintf(tempOut, "%s/replay_pass4_%s.root",getenv("REPLAY_PASS4"), argv[1]);
  TFile *fileOut = new TFile(tempOut,"RECREATE");
  TTree *Tout = new TTree("pass4", "pass4");

  // Variables
  Tout->Branch("pmt0_pass4", &pmt_pass4[0], "pmt0_pass4/D");
  Tout->Branch("pmt1_pass4", &pmt_pass4[1], "pmt1_pass4/D");
  Tout->Branch("pmt2_pass4", &pmt_pass4[2], "pmt2_pass4/D");
  Tout->Branch("pmt3_pass4", &pmt_pass4[3], "pmt3_pass4/D");
  Tout->Branch("pmt4_pass4", &pmt_pass4[4], "pmt4_pass4/D");
  Tout->Branch("pmt5_pass4", &pmt_pass4[5], "pmt5_pass4/D");
  Tout->Branch("pmt6_pass4", &pmt_pass4[6], "pmt6_pass4/D");
  Tout->Branch("pmt7_pass4", &pmt_pass4[7], "pmt7_pass4/D");

  Tout->Branch("AnodeE", &AnodeE, "AnodeE/D");
  Tout->Branch("AnodeW", &AnodeW, "AnodeW/D");

  Tout->Branch("timeE", &timeE, "timeE/D");
  Tout->Branch("timeW", &timeW, "timeW/D");
  Tout->Branch("timeE_BB", &timeE_BB, "timeE_BB/D");
  Tout->Branch("timeW_BB", &timeW_BB, "timeW_BB/D");
  Tout->Branch("UBtime", &UBtime, "UBtime/D");
  Tout->Branch("UBtime_BB", &UBtime_BB, "UBtime_BB/D");
  Tout->Branch("twoFoldE", &twoFoldE, "twoFoldE/D");
  Tout->Branch("twoFoldW", &twoFoldW, "twoFoldW/D");

  //Branch with leaves to store the Evis and weights for all 8 PMTs
  Tout->Branch("pmt_Evis", &pmt_Evis,
	       "Evis0/D:Evis1:Evis2:Evis3:Evis4:Evis5:Evis6:Evis7:weight0:weight1:weight2:weight3:weight4:weight5:weight6:weight7");

  Tout->Branch("xE_pass4", &xE_pass4, "xE_pass4/D");
  Tout->Branch("yE_pass4", &yE_pass4, "yE_pass4/D");
  Tout->Branch("xW_pass4", &xW_pass4, "xW_pass4/D");
  Tout->Branch("yW_pass4", &yW_pass4, "yW_pass4/D");
  
  Tout->Branch("EreconTot", &EreconTot, "EreconTot/D");
  Tout->Branch("EreconE", &EreconE, "EreconE/D");
  Tout->Branch("EreconW", &EreconW, "EreconW/D");

  Tout->Branch("PID_pass4",  &PID_pass4,  "PID_pass4/I");
  Tout->Branch("type_pass4", &type_pass4, "type_pass4/I");
  Tout->Branch("side_pass4", &side_pass4, "side_pass4/I");
  Tout->Branch("posError_pass4", &posError_pass4, "posError_pass4/I");

  // Open input ntuple from pass3
  char tempIn[500];
  sprintf(tempIn, "%s/replay_pass3_%s.root", getenv("REPLAY_PASS3"),argv[1]);

  TFile *fileIn = new TFile(tempIn, "READ");
  TTree *Tin = (TTree*)(fileIn->Get("pass3"));

  // Variables
  Tin->SetBranchAddress("pmt0_pass3", &pmt_pass3[0]);
  Tin->SetBranchAddress("pmt1_pass3", &pmt_pass3[1]);
  Tin->SetBranchAddress("pmt2_pass3", &pmt_pass3[2]);
  Tin->SetBranchAddress("pmt3_pass3", &pmt_pass3[3]);
  Tin->SetBranchAddress("pmt4_pass3", &pmt_pass3[4]);
  Tin->SetBranchAddress("pmt5_pass3", &pmt_pass3[5]);
  Tin->SetBranchAddress("pmt6_pass3", &pmt_pass3[6]);
  Tin->SetBranchAddress("pmt7_pass3", &pmt_pass3[7]);
  
  Tin->SetBranchAddress("AnodeE", &AnodeE);
  Tin->SetBranchAddress("AnodeW", &AnodeW);

  Tin->SetBranchAddress("timeE", &timeE);
  Tin->SetBranchAddress("timeW", &timeW);
  Tin->SetBranchAddress("timeE_BB", &timeE_BB);
  Tin->SetBranchAddress("timeW_BB", &timeW_BB);
  Tin->SetBranchAddress("UBtime", &UBtime);
  Tin->SetBranchAddress("UBtime_BB", &UBtime_BB);
  Tin->SetBranchAddress("twoFoldE", &twoFoldE);
  Tin->SetBranchAddress("twoFoldW", &twoFoldW);

  Tin->SetBranchAddress("xE_pass3", &xE_pass3);
  Tin->SetBranchAddress("yE_pass3", &yE_pass3);
  Tin->SetBranchAddress("xW_pass3", &xW_pass3);
  Tin->SetBranchAddress("yW_pass3", &yW_pass3);

  Tin->SetBranchAddress("PID_pass3",  &PID_pass3);
  Tin->SetBranchAddress("type_pass3", &type_pass3);
  Tin->SetBranchAddress("side_pass3", &side_pass3);
  Tin->SetBranchAddress("posError_pass3", &posError_pass3);

  //Open input ntuple from pass2 for the non-position corrected ADC value for calculating the weight
  sprintf(tempIn, "%s/replay_pass2_%s.root", getenv("REPLAY_PASS2"),argv[1]);

  TFile *fileIn2 = new TFile(tempIn, "READ");
  TTree *Tin2 = (TTree*)(fileIn2->Get("pass2"));

  // Variables
  Tin2->SetBranchAddress("pmt0_pass2", &pmt_pass2[0]);
  Tin2->SetBranchAddress("pmt1_pass2", &pmt_pass2[1]);
  Tin2->SetBranchAddress("pmt2_pass2", &pmt_pass2[2]);
  Tin2->SetBranchAddress("pmt3_pass2", &pmt_pass2[3]);
  Tin2->SetBranchAddress("pmt4_pass2", &pmt_pass2[4]);
  Tin2->SetBranchAddress("pmt5_pass2", &pmt_pass2[5]);
  Tin2->SetBranchAddress("pmt6_pass2", &pmt_pass2[6]);
  Tin2->SetBranchAddress("pmt7_pass2", &pmt_pass2[7]);

  Tin2->SetBranchAddress("xE_pass2", &xE_pass2);
  Tin2->SetBranchAddress("yE_pass2", &yE_pass2);
  Tin2->SetBranchAddress("xW_pass2", &xW_pass2);
  Tin2->SetBranchAddress("yW_pass2", &yW_pass2);

  Tin2->SetBranchAddress("PID_pass2",  &PID_pass2);
  Tin2->SetBranchAddress("type_pass2", &type_pass2);
  Tin2->SetBranchAddress("side_pass2", &side_pass2);
  Tin2->SetBranchAddress("posError_pass2", &posError_pass2);

  int nEvents = Tin->GetEntriesFast();
  cout << "... Processing nEvents = " << nEvents << endl;

  // Loop over events
  for (int i=0; i<nEvents; i++) {
    Tin->GetEvent(i);
    Tin2->GetEvent(i);

    //Calculate Evis for each event in each PMT using the linearity curve determined in calibration
    pmt_Evis.Evis0 = linearityCurve[0][0] + linearityCurve[0][1]*pmt_pass3[0] + linearityCurve[0][2]*pmt_pass3[0]*pmt_pass3[0];
    pmt_Evis.Evis1 = linearityCurve[1][0] + linearityCurve[1][1]*pmt_pass3[1] + linearityCurve[1][2]*pmt_pass3[1]*pmt_pass3[1];
    pmt_Evis.Evis2 = linearityCurve[2][0] + linearityCurve[2][1]*pmt_pass3[2] + linearityCurve[2][2]*pmt_pass3[2]*pmt_pass3[2];
    pmt_Evis.Evis3 = linearityCurve[3][0] + linearityCurve[3][1]*pmt_pass3[3] + linearityCurve[3][2]*pmt_pass3[3]*pmt_pass3[3];
    pmt_Evis.Evis4 = linearityCurve[4][0] + linearityCurve[4][1]*pmt_pass3[4] + linearityCurve[4][2]*pmt_pass3[4]*pmt_pass3[4];
    pmt_Evis.Evis5 = linearityCurve[5][0] + linearityCurve[5][1]*pmt_pass3[5] + linearityCurve[5][2]*pmt_pass3[5]*pmt_pass3[5];
    pmt_Evis.Evis6 = linearityCurve[6][0] + linearityCurve[6][1]*pmt_pass3[6] + linearityCurve[6][2]*pmt_pass3[6]*pmt_pass3[6];
    pmt_Evis.Evis7 = linearityCurve[7][0] + linearityCurve[7][1]*pmt_pass3[7] + linearityCurve[7][2]*pmt_pass3[7]*pmt_pass3[7];

    //Now map each Evis value to a true value using EQ2Etrue relationship as was determined in simulation
    Etrue[0] = EQ2Etrue[0][0]+EQ2Etrue[0][1]*(pmt_Evis.Evis0)+EQ2Etrue[0][2]*(pmt_Evis.Evis0)*(pmt_Evis.Evis0);
    Etrue[1] = EQ2Etrue[1][0]+EQ2Etrue[1][1]*(pmt_Evis.Evis1)+EQ2Etrue[1][2]*(pmt_Evis.Evis1)*(pmt_Evis.Evis1);
    Etrue[2] = EQ2Etrue[2][0]+EQ2Etrue[2][1]*(pmt_Evis.Evis2)+EQ2Etrue[2][2]*(pmt_Evis.Evis2)*(pmt_Evis.Evis2);
    Etrue[3] = EQ2Etrue[3][0]+EQ2Etrue[3][1]*(pmt_Evis.Evis3)+EQ2Etrue[3][2]*(pmt_Evis.Evis3)*(pmt_Evis.Evis3);
    Etrue[4] = EQ2Etrue[4][0]+EQ2Etrue[4][1]*(pmt_Evis.Evis4)+EQ2Etrue[4][2]*(pmt_Evis.Evis4)*(pmt_Evis.Evis4);
    Etrue[5] = EQ2Etrue[5][0]+EQ2Etrue[5][1]*(pmt_Evis.Evis5)+EQ2Etrue[5][2]*(pmt_Evis.Evis5)*(pmt_Evis.Evis5);
    Etrue[6] = EQ2Etrue[6][0]+EQ2Etrue[6][1]*(pmt_Evis.Evis6)+EQ2Etrue[6][2]*(pmt_Evis.Evis6)*(pmt_Evis.Evis6);
    Etrue[7] = EQ2Etrue[7][0]+EQ2Etrue[7][1]*(pmt_Evis.Evis7)+EQ2Etrue[7][2]*(pmt_Evis.Evis7)*(pmt_Evis.Evis7);

    if (pmtQuality[0] && pmt_Evis.Evis0>0. && (side_pass3==0 || type_pass3==1)) {
      Double_t N = pmt_pass2[0]*nPE_per_channel[0];
      Double_t f = sqrt(N)/N;
      pmt_Evis.weight0 = 1/(pmt_Evis.Evis0*pmt_Evis.Evis0*f*f);
    }
    else {pmt_Evis.weight0=0.;}

    if (pmtQuality[1] && pmt_Evis.Evis1>0. && (side_pass3==0 || type_pass3==1)) {
      Double_t N = pmt_pass2[1]*nPE_per_channel[1];
      Double_t f = sqrt(N)/N;
      pmt_Evis.weight1 = 1/(pmt_Evis.Evis1*pmt_Evis.Evis1*f*f);
    }
    else {pmt_Evis.weight1=0.;}

    if (pmtQuality[2] && pmt_Evis.Evis2>0. && (side_pass3==0 || type_pass3==1)) {
      Double_t N = pmt_pass2[2]*nPE_per_channel[2];
      Double_t f = sqrt(N)/N;
      pmt_Evis.weight2 = 1/(pmt_Evis.Evis2*pmt_Evis.Evis2*f*f);
    }
    else {pmt_Evis.weight2=0.;}

    if (pmtQuality[3] && pmt_Evis.Evis3>0. && (side_pass3==0 || type_pass3==1)) {
      Double_t N = pmt_pass2[3]*nPE_per_channel[3];
      Double_t f = sqrt(N)/N;
      pmt_Evis.weight3 = 1/(pmt_Evis.Evis3*pmt_Evis.Evis3*f*f);
    }
    else {pmt_Evis.weight3=0.;}

    if (pmtQuality[4] && pmt_Evis.Evis4>0. && (side_pass3==1 || type_pass3==1)) {
      Double_t N = pmt_pass2[4]*nPE_per_channel[4];
      Double_t f = sqrt(N)/N;
      pmt_Evis.weight4 = 1/(pmt_Evis.Evis4*pmt_Evis.Evis4*f*f);
    }
    else {pmt_Evis.weight4=0.;}

    if (pmtQuality[5] && pmt_Evis.Evis5>0. && (side_pass3==1 || type_pass3==1)) {
      Double_t N = pmt_pass2[5]*nPE_per_channel[5];
      Double_t f = sqrt(N)/N;
      pmt_Evis.weight5 = 1/(pmt_Evis.Evis5*pmt_Evis.Evis5*f*f);
    }
    else {pmt_Evis.weight5=0.;}

    if (pmtQuality[6] && pmt_Evis.Evis6>0. && (side_pass3==1 || type_pass3==1)) {
      Double_t N = pmt_pass2[6]*nPE_per_channel[6];
      Double_t f = sqrt(N)/N;
      pmt_Evis.weight6 = 1/(pmt_Evis.Evis6*pmt_Evis.Evis6*f*f);
    }
    else {pmt_Evis.weight6=0.;}

    if (pmtQuality[7] && pmt_Evis.Evis7>0. && (side_pass3==1 || type_pass3==1)) {
      Double_t N = pmt_pass2[7]*nPE_per_channel[7];
      Double_t f = sqrt(N)/N;
      pmt_Evis.weight7 = 1/(pmt_Evis.Evis7*pmt_Evis.Evis7*f*f);
    }
    else {pmt_Evis.weight7=0.;}

    //East side EvisE
    if (side_pass3==0 || type_pass3==1) {
      EreconE = (pmt_Evis.weight0*Etrue[0]+pmt_Evis.weight1*Etrue[1]+pmt_Evis.weight2*Etrue[2]+pmt_Evis.weight3*Etrue[3])/(pmt_Evis.weight0+pmt_Evis.weight1+pmt_Evis.weight2+pmt_Evis.weight3);}
    else EreconE=0.;
    //West Side EvisW
    if (side_pass3==1 || type_pass3==1) {
      EreconW = (pmt_Evis.weight4*Etrue[4]+pmt_Evis.weight5*Etrue[5]+pmt_Evis.weight6*Etrue[6]+pmt_Evis.weight7*Etrue[7])/(pmt_Evis.weight4+pmt_Evis.weight5+pmt_Evis.weight6+pmt_Evis.weight7);}
    else EreconW=0.;
    
    EreconTot = EreconE+EreconW;

    if (i<20) {
      cout << Etrue[4] << " " << Etrue[5] << " " << Etrue[6] << " " << Etrue[7] << endl;
      cout << pmt_Evis.weight4 << " " << pmt_Evis.weight5 << " " << pmt_Evis.weight6 << " " << pmt_Evis.weight7 << endl;
    }
    
    // Pass other variables from pass3 to pass4
    pmt_pass4[0] = pmt_pass3[0];
    pmt_pass4[1] = pmt_pass3[1];
    pmt_pass4[2] = pmt_pass3[2];
    pmt_pass4[3] = pmt_pass3[3];
    pmt_pass4[4] = pmt_pass3[4];
    pmt_pass4[5] = pmt_pass3[5];
    pmt_pass4[6] = pmt_pass3[6];
    pmt_pass4[7] = pmt_pass3[7];

    xE_pass4 = xE_pass3;
    yE_pass4 = yE_pass3;
    xW_pass4 = xW_pass3;
    yW_pass4 = yW_pass3;

    PID_pass4 = PID_pass3;
    type_pass4 = type_pass3;
    side_pass4 = side_pass3;
    posError_pass4 = posError_pass3;

    Tout->Fill();
  }

  // Write output ntuple
  fileOut->Write();
  fileOut->Close();

  return 0;
}
