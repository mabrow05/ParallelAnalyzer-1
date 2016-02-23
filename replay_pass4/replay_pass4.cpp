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
#include "DataTree.hh"
#include "runInfo.h"

using namespace std;

int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);  

  Int_t nPMTs = 8;
  Int_t nParams = 6; //cubic fit plus two terms for the linear piece wise term where p0->p3 are for cubic and 
                     // and p4 is the where the low energy fit becomes linear and p5 is the slope of this line
  
  // Run number integer
  cout << "Run " << argv[1] << " ..." << endl;
  cout << "... Applying Energy Calibration ..." << endl;
  //istringstream ss(argv[1]);
  int runNumber = atoi(argv[1]);
  //ss >> runNumber;
  unsigned int calibrationPeriod = getSrcRunPeriod(runNumber);
  // Determine linearity curve to use
  char tempFileLinearityCurve[500];
  
  sprintf(tempFileLinearityCurve, "../linearity_curves/lin_curves_srcCal_Period_%i.dat",calibrationPeriod);
  cout << "... Reading: " << tempFileLinearityCurve << endl;

   //Setup array to hold linearity curves
  Double_t linearityCurve[nPMTs][nParams];

  // Read linearity curve
  cout << "Reading in linearity curve:\n";
  cout << "p0\tp1\tp2\tp3\tp4\tp5\n";
  ifstream fileLinearityCurve(tempFileLinearityCurve);
  Double_t p[6];//p0,p1,p2,p3,p4,p5;
  Int_t i=0;
  while (fileLinearityCurve >> p[0] >> p[1] >> p[2] >> p[3] >> p[4] >> p[5]) {
    Int_t ii=0;
    while(ii<nParams) {linearityCurve[i][ii] = p[ii]; ii++;}
    i++;
    cout << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " " << p[4] << " " << p[5] << endl;
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
  //sprintf(tempOut, "replay_pass4_%s.root", argv[1]);
  DataTree *t = new DataTree();
  t->makeOutputTree(std::string(tempOut),"pass4");

  char tempIn[500];
  sprintf(tempIn, "%s/replay_pass3_%s.root", getenv("REPLAY_PASS3"),argv[1]);
  //sprintf(tempIn, "../replay_pass3/replay_pass3_%s.root",argv[1]);
  t->setupInputTree(std::string(tempIn),"pass3");
  
  
  int nEvents = t->getEntries();
  cout << "... Processing nEvents = " << nEvents << endl;

  int pp=0; // counter for events to print to screen

  // Loop over events
  for (int i=0; i<nEvents; i++) {
    t->getEvent(i);
    
    //Calculate Evis for each event in each PMT using the linearity curve determined in calibration
    //if (t->ScintE.q1>linearityCurve[0][4]) {
      t->ScintE.e1 = linearityCurve[0][0] + linearityCurve[0][1]*t->ScintE.q1 
	+ linearityCurve[0][2]*t->ScintE.q1*t->ScintE.q1 + linearityCurve[0][3]*t->ScintE.q1*t->ScintE.q1*t->ScintE.q1;
      //}
      //else if (t->ScintE.q1>0.) t->ScintE.e1 = linearityCurve[0][5]*t->ScintE.q1;
      //else t->ScintE.e1=0.;

      //if (t->ScintE.q2>linearityCurve[1][4]) {
      t->ScintE.e2 = linearityCurve[1][0] + linearityCurve[1][1]*t->ScintE.q2 
	+ linearityCurve[1][2]*t->ScintE.q2*t->ScintE.q2 + linearityCurve[1][3]*t->ScintE.q2*t->ScintE.q2*t->ScintE.q2;
      //}
      //else if (t->ScintE.q2>0.) t->ScintE.e2 = linearityCurve[1][5]*t->ScintE.q2;
      //else t->ScintE.e2=0.;

      //if (t->ScintE.q3>linearityCurve[2][4]) {
      t->ScintE.e3 = linearityCurve[2][0] + linearityCurve[2][1]*t->ScintE.q3 
	+ linearityCurve[2][2]*t->ScintE.q3*t->ScintE.q3 + linearityCurve[2][3]*t->ScintE.q3*t->ScintE.q3*t->ScintE.q3;
      //}
    //else if (t->ScintE.q3>0.) t->ScintE.e3 = linearityCurve[2][5]*t->ScintE.q3;
    //else t->ScintE.e3=0.;

    //if (t->ScintE.q4>linearityCurve[3][4]) {
      t->ScintE.e4 = linearityCurve[3][0] + linearityCurve[3][1]*t->ScintE.q4 
	+ linearityCurve[3][2]*t->ScintE.q4*t->ScintE.q4 + linearityCurve[3][3]*t->ScintE.q4*t->ScintE.q4*t->ScintE.q4;
      //}
    //else if (t->ScintE.q4>0.) t->ScintE.e4 = linearityCurve[3][5]*t->ScintE.q4;
    //else t->ScintE.e4=0.;

    //if (t->ScintW.q1>linearityCurve[4][4]) {
      t->ScintW.e1 = linearityCurve[4][0] + linearityCurve[4][1]*t->ScintW.q1 
	+ linearityCurve[4][2]*t->ScintW.q1*t->ScintW.q1 + linearityCurve[4][3]*t->ScintW.q1*t->ScintW.q1*t->ScintW.q1;
      //}
    //else if (t->ScintW.q1>0.) t->ScintW.e1 = linearityCurve[4][5]*t->ScintW.q1;
    //else t->ScintW.e1=0.;

    //if (t->ScintW.q2>linearityCurve[5][4]) {
      t->ScintW.e2 = linearityCurve[5][0] + linearityCurve[5][1]*t->ScintW.q2 
	+ linearityCurve[5][2]*t->ScintW.q2*t->ScintW.q2 + linearityCurve[5][3]*t->ScintW.q2*t->ScintW.q2*t->ScintW.q2;
      //}
    //else if (t->ScintW.q2>0.) t->ScintW.e2 = linearityCurve[5][5]*t->ScintW.q2;
    //else t->ScintW.e2=0.;

    //if (t->ScintW.q3>linearityCurve[6][4]) {
      t->ScintW.e3 = linearityCurve[6][0] + linearityCurve[6][1]*t->ScintW.q3 
	+ linearityCurve[6][2]*t->ScintW.q3*t->ScintW.q3 + linearityCurve[6][3]*t->ScintW.q3*t->ScintW.q3*t->ScintW.q3;
      //}
    //else if (t->ScintW.q3>0.) t->ScintW.e3 = linearityCurve[6][5]*t->ScintW.q3;
    //else t->ScintW.e3=0.;

    //if (t->ScintW.q4>linearityCurve[7][4]) {
      t->ScintW.e4 = linearityCurve[7][0] + linearityCurve[7][1]*t->ScintW.q4 
	+ linearityCurve[7][2]*t->ScintW.q4*t->ScintW.q4 + linearityCurve[7][3]*t->ScintW.q4*t->ScintW.q4*t->ScintW.q4;
      //}
    //else if (t->ScintW.q4>0.) t->ScintW.e4 = linearityCurve[7][5]*t->ScintW.q4;
    //else t->ScintW.e4=0.;

    //Now map each Evis value to a true value using EQ2Etrue relationship as was determined in simulation
    /*Etrue[0] = EQ2Etrue[0][0]+EQ2Etrue[0][1]*(t->ScintE.e1)+EQ2Etrue[0][2]*(t->ScintE.e1)*(t->ScintE.e1);
    Etrue[1] = EQ2Etrue[1][0]+EQ2Etrue[1][1]*(t->ScintE.e2)+EQ2Etrue[1][2]*(t->ScintE.e2)*(t->ScintE.e2);
    Etrue[2] = EQ2Etrue[2][0]+EQ2Etrue[2][1]*(t->ScintE.e3)+EQ2Etrue[2][2]*(t->ScintE.e3)*(t->ScintE.e3);
    Etrue[3] = EQ2Etrue[3][0]+EQ2Etrue[3][1]*(t->ScintE.e4)+EQ2Etrue[3][2]*(t->ScintE.e4)*(t->ScintE.e4);
    Etrue[4] = EQ2Etrue[4][0]+EQ2Etrue[4][1]*(t->ScintW.e1)+EQ2Etrue[4][2]*(t->ScintW.e1)*(t->ScintW.e1);
    Etrue[5] = EQ2Etrue[5][0]+EQ2Etrue[5][1]*(t->ScintW.e2)+EQ2Etrue[5][2]*(t->ScintW.e2)*(t->ScintW.e2);
    Etrue[6] = EQ2Etrue[6][0]+EQ2Etrue[6][1]*(t->ScintW.e3)+EQ2Etrue[6][2]*(t->ScintW.e3)*(t->ScintW.e3);
    Etrue[7] = EQ2Etrue[7][0]+EQ2Etrue[7][1]*(t->ScintW.e4)+EQ2Etrue[7][2]*(t->ScintW.e4)*(t->ScintW.e4);
    */
    Double_t lowestPMTenergy = 1.; //This is for setting an upper limit on the weight since
                                        // at E=0, the weight goes to infiniti
    Double_t lowestADC = 0.001;//20.;

    if (pmtQuality[0]) { // && t->ScintE.e1>0.) { 
      if (t->ScintE.e1>lowestPMTenergy && t->ScintE.q1>lowestADC) {
	t->ScintE.nPE1 = t->ScintE.q1*nPE_per_channel[0];
	t->ScintE.de1 = t->ScintE.e1/sqrt(t->ScintE.nPE1);
	pmt_Evis.weight0 = t->ScintE.nPE1/(t->ScintE.e1*t->ScintE.e1);
      }
      else pmt_Evis.weight0 = lowestADC*nPE_per_channel[0]/(lowestPMTenergy*lowestPMTenergy);
    }
    else {pmt_Evis.weight0=0.;}

    if (pmtQuality[1]) { // && t->ScintE.e2>0.) { 
      if (t->ScintE.e2>lowestPMTenergy && t->ScintE.q2>lowestADC) {
	t->ScintE.nPE2 = t->ScintE.q2*nPE_per_channel[1];
	t->ScintE.de2 = t->ScintE.e2/sqrt(t->ScintE.nPE2);
	pmt_Evis.weight1 = t->ScintE.nPE2/(t->ScintE.e2*t->ScintE.e2);
      }
      else pmt_Evis.weight1 = lowestADC*nPE_per_channel[1]/(lowestPMTenergy*lowestPMTenergy);
    }
    else {pmt_Evis.weight1=0.;}

    if (pmtQuality[2]) { // && t->ScintE.e3>0.) { 
      if (t->ScintE.e3>lowestPMTenergy && t->ScintE.q3>lowestADC) {
	t->ScintE.nPE3 = t->ScintE.q3*nPE_per_channel[2];
	t->ScintE.de3 = t->ScintE.e3/sqrt(t->ScintE.nPE3);
	pmt_Evis.weight2 = t->ScintE.nPE3/(t->ScintE.e3*t->ScintE.e3);
      }
      else pmt_Evis.weight2 = lowestADC*nPE_per_channel[2]/(lowestPMTenergy*lowestPMTenergy);
    }
    else {pmt_Evis.weight2=0.;}

    if (pmtQuality[3]) { // && t->ScintE.e4>0.) { 
      if (t->ScintE.e4>lowestPMTenergy && t->ScintE.q4>lowestADC) {
	t->ScintE.nPE4 = t->ScintE.q4*nPE_per_channel[3];
	t->ScintE.de4 = t->ScintE.e4/sqrt(t->ScintE.nPE4);
	pmt_Evis.weight3 = t->ScintE.nPE4/(t->ScintE.e4*t->ScintE.e4);
      }
      else pmt_Evis.weight3 = lowestADC*nPE_per_channel[3]/(lowestPMTenergy*lowestPMTenergy);
    }
    else {pmt_Evis.weight3=0.;}

    if (pmtQuality[4]) { // && t->ScintW.e1>0.) { 
      if (t->ScintW.e1>lowestPMTenergy && t->ScintW.q1>lowestADC) {
	t->ScintW.nPE1 = t->ScintW.q1*nPE_per_channel[4];
	t->ScintW.de1 = t->ScintW.e1/sqrt(t->ScintW.nPE1);
	pmt_Evis.weight4 = t->ScintW.nPE1/(t->ScintW.e1*t->ScintW.e1);
      }
      else pmt_Evis.weight4 = lowestADC*nPE_per_channel[4]/(lowestPMTenergy*lowestPMTenergy);
    }
    else {pmt_Evis.weight4=0.;}

    if (pmtQuality[5]) { // && t->ScintW.e2>0.) {
      if (t->ScintW.e2>lowestPMTenergy && t->ScintW.q2>lowestADC) {
	t->ScintW.nPE2 = t->ScintW.q2*nPE_per_channel[5];
	t->ScintW.de2 = t->ScintW.e2/sqrt(t->ScintW.nPE2);
	pmt_Evis.weight5 = t->ScintW.nPE2/(t->ScintW.e2*t->ScintW.e2);
      }
      else pmt_Evis.weight5 = lowestADC*nPE_per_channel[5]/(lowestPMTenergy*lowestPMTenergy);
    }
    else {pmt_Evis.weight5=0.;}

    if (pmtQuality[6]) { // && t->ScintW.e3>0.) {
      if (t->ScintW.e3>lowestPMTenergy && t->ScintW.q3>lowestADC) {
	t->ScintW.nPE3 = t->ScintW.q3*nPE_per_channel[6];
	t->ScintW.de3 = t->ScintW.e3/sqrt(t->ScintW.nPE3);
	pmt_Evis.weight6 = t->ScintW.nPE3/(t->ScintW.e3*t->ScintW.e3);
      }
      else pmt_Evis.weight6 = lowestADC*nPE_per_channel[6]/(lowestPMTenergy*lowestPMTenergy);
    }
    else {pmt_Evis.weight6=0.;}

    if (pmtQuality[7]) { // && t->ScintW.e4>0.) {
      if (t->ScintW.e4>lowestPMTenergy && t->ScintW.q4>lowestADC) {
	t->ScintW.nPE4 = t->ScintW.q4*nPE_per_channel[7];
	t->ScintW.de4 = t->ScintW.e4/sqrt(t->ScintW.nPE4);
	pmt_Evis.weight7 = t->ScintW.nPE4/(t->ScintW.e4*t->ScintW.e4);
      }
      else pmt_Evis.weight7 = lowestADC*nPE_per_channel[7]/(lowestPMTenergy*lowestPMTenergy);
    }
    else {pmt_Evis.weight7=0.;}

    //East side EvisE
    //if (!(pmt_Evis.weight0==0. && pmt_Evis.weight1==0. && pmt_Evis.weight2==0. && pmt_Evis.weight3==0.)) {
      t->ScintE.energy = t->EvisE = (pmt_Evis.weight0*t->ScintE.e1+pmt_Evis.weight1*t->ScintE.e2+pmt_Evis.weight2*t->ScintE.e3+pmt_Evis.weight3*t->ScintE.e4)/(pmt_Evis.weight0+pmt_Evis.weight1+pmt_Evis.weight2+pmt_Evis.weight3);
      t->ScintE.denergy = 1./sqrt(pmt_Evis.weight0+pmt_Evis.weight1+pmt_Evis.weight2+pmt_Evis.weight3);
      //}
      //else {
      //t->ScintE.energy = t->EvisE = 0.;
      //t->ScintE.denergy = 0.;
      // }
//else EvisE=0.;
    //West Side EvisW
    //if (!(pmt_Evis.weight4==0. && pmt_Evis.weight5==0. && pmt_Evis.weight6==0. && pmt_Evis.weight7==0.)) {
      t->ScintW.energy = t->EvisW = (pmt_Evis.weight4*t->ScintW.e1+pmt_Evis.weight5*t->ScintW.e2+pmt_Evis.weight6*t->ScintW.e3+pmt_Evis.weight7*t->ScintW.e4)/(pmt_Evis.weight4+pmt_Evis.weight5+pmt_Evis.weight6+pmt_Evis.weight7);
      t->ScintW.denergy = 1./sqrt(pmt_Evis.weight4+pmt_Evis.weight5+pmt_Evis.weight6+pmt_Evis.weight7);
      //}
      //else {
      //t->ScintW.energy = t->EvisW = 0.;
      //t->ScintW.denergy = 0.;
      // }
      //else EvisW=0.;
    
      t->Erecon = (t->EvisW>0. && t->EvisE>0.) ? t->EvisE+t->EvisW : (t->EvisW>0. ? t->EvisW : (t->EvisE>0. ? t->EvisE : 0.));

    
    //if (t->Side==1) cout << pmt_Evis.weight4 << " " << pmt_Evis.weight5 << " " << pmt_Evis.weight6 << " " << pmt_Evis.weight7 << endl;
    //if (t->ScintW.energy>150.) {
    // while (pp<20000) {
    //	cout << pmt_Evis.weight4 << " " << pmt_Evis.weight5 << " " << pmt_Evis.weight6 << " " << pmt_Evis.weight7 << endl;
    //	pp++;
    //}
    //}
    
    if (t->Erecon>0.) t->fillOutputTree();
  }

  // Write output ntuple
  t->writeOutputFile();
  delete t;
  

  return 0;
}


/////////// OLD INPUT OUTPUT METHOD /////////////////
/*TFile *fileOut = new TFile(tempOut,"RECREATE");
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

  int xeRC,yeRC,xwRC,ywRC;	
  Tout->Branch("xeRC", &xeRC, "xeRC/I"); //x east response class. 
  Tout->Branch("yeRC", &yeRC, "yeRC/I"); //y east response class... 
  Tout->Branch("xwRC", &xwRC, "xwRC/I");
  Tout->Branch("ywRC", &ywRC, "ywRC/I");

  Tout->Branch("EvisTot", &EvisTot, "EvisTot/D");
  Tout->Branch("EvisE", &EvisE, "EvisE/D");
  Tout->Branch("EvisW", &EvisW, "EvisW/D");
  
  //Tout->Branch("EreconTot", &EreconTot, "EreconTot/D");
  //Tout->Branch("EreconE", &EreconE, "EreconE/D");
  //Tout->Branch("EreconW", &EreconW, "EreconW/D");

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

  //Adding in wirechamber class variable
  Tin->SetBranchAddress("xeRC",&xeRC);
  Tin->SetBranchAddress("yeRC",&yeRC);
  Tin->SetBranchAddress("xwRC",&xwRC);
  Tin->SetBranchAddress("ywRC",&ywRC);

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


*/
