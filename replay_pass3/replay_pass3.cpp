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
#include "runInfo.h"
#include "DataTree.hh"
#include "posMapReader.h"
#include "positionMapHandler.hh"
#include "calibrationTools.hh"

#include "replay_pass2.h"
#include "replay_pass3.h"

std::vector < std::vector < std::vector <double> > > getEQ2EtrueParams(int runNumber) {
  ifstream infile;
  Char_t temp[200];
  if (runNumber<20000) {
    sprintf(temp,"%s/simulation_comparison/EQ2EtrueConversion/2011-2012_EQ2EtrueFitParams.dat", getenv("ANALYSIS_CODE"));
    infile.open(temp);
  }
  else {
    sprintf(temp,"%s/simulation_comparison/EQ2EtrueConversion/2012-2013_EQ2EtrueFitParams.dat", getenv("ANALYSIS_CODE"));
    infile.open(temp);
  }
  
  std::vector < std::vector < std::vector < double > > > params;
  params.resize(2,std::vector < std::vector < double > > (3, std::vector < double > (6,0.)));

  char holdType[10];
  int side=0, type=0;
  while (infile >> holdType >> params[side][type][0] >> params[side][type][1] >> params[side][type][2] >> params[side][type][3] >> params[side][type][4] >> params[side][type][5]) { 
    cout << holdType << " " << params[side][type][0] << " " << params[side][type][1] << " " << params[side][type][2] << " " << params[side][type][3] << " " << params[side][type][4] << " " << params[side][type][5] << endl;
    type+=1;
    if (type==3) {type=0; side=1;}
  }
  return params;
};

vector < Double_t > GetAlphaValues(Int_t runPeriod)
{
  Char_t temp[500];
  vector < Double_t > alphas (8,0.);
  sprintf(temp,"%s/simulation_comparison/nPE_per_keV/nPE_per_keV_%i.dat",getenv("ANALYSIS_CODE"),runPeriod);
  ifstream infile;
  infile.open(temp);
  Int_t i = 0;

  while (infile >> alphas[i]) { std::cout << alphas[i] << std::endl; i++; }
  return alphas;
};

vector <Int_t> getPMTQuality(Int_t runNumber) {
  //Read in PMT quality file
  cout << "Reading in PMT Quality file ...\n";
  vector <Int_t>  pmtQuality (8,0);
  Char_t temp[200];
  sprintf(temp,"%s/residuals/PMT_runQuality_master.dat",getenv("ANALYSIS_CODE")); 
  ifstream pmt;
  std::cout << temp << std::endl;
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
  return pmtQuality;
};


using namespace std;

int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);
  
  int runNumber = atoi(argv[1]);
  int nPMT = 8;
  int nParams = 3; //takes a quadratic (but quadratic term set to zero)

  // Run number integer
  cout << "Run " << runNumber << " ..." << endl;
  cout << "... Applying Calibration ..." << endl;

  unsigned int calibrationPeriod = getSrcRunPeriod(runNumber);

  //Get the linearity curve
  LinearityCurve linearityCurve(calibrationPeriod);

  //Load the simulated relationship between EQ and Etrue
  vector < vector < vector < double > > > EQ2Etrue = getEQ2EtrueParams(runNumber);
 
  //Read in PMT quality file
  std::vector <Int_t> pmtQuality = getPMTQuality(runNumber);

  //Get values for nPE/keV...
  std::vector <Double_t> alpha = GetAlphaValues(calibrationPeriod);
  
  //Reading position map...
  UInt_t XePeriod = getXeRunPeriod(runNumber); // Get the proper Xe run period for the Trigger functions
  //GetPositionMap(XePeriod);
  PositionMap posmap(5.0);
  posmap.readPositionMap(XePeriod);

  // DataTree structure
  DataTree *t = new DataTree();

  // Open output ntuple
  char tempOut[500];
  sprintf(tempOut, "%s/replay_pass3_%s.root",getenv("REPLAY_PASS3"), argv[1]);
  //sprintf(tempOut, "replay_pass3_%s.root", argv[1]);
  t->makeOutputTree(std::string(tempOut),"pass3");

  // Input ntuple
  char tempIn[500];
  sprintf(tempIn, "%s/replay_pass2_%s.root", getenv("REPLAY_PASS2"),argv[1]);
  //sprintf(tempIn, "../replay_pass2/replay_pass2_%s.root",argv[1]);
  t->setupInputTree(std::string(tempIn),"pass2");
 
  int nEvents = t->getEntries();
  cout << "... Processing nEvents = " << nEvents << endl;

  vector < vector <Int_t> > gridPoint;
  vector < Double_t > eta;

  // Loop over events
  for (int i=0; i<nEvents; i++) {
    t->getEvent(i);

    eta = posmap.getInterpolatedEta(t->xE.center, t->yE.center, t->xW.center, t->yW.center);

    
    t->ScintE.e1 = eta[0]>0. ? linearityCurve.applyLinCurve(0,t->ScintE.q1) / eta[0] : 0.;
    t->ScintE.e2 = eta[1]>0. ? linearityCurve.applyLinCurve(1,t->ScintE.q2) / eta[1] : 0.;
    t->ScintE.e3 = eta[2]>0. ? linearityCurve.applyLinCurve(2,t->ScintE.q3) / eta[2] : 0.;
    t->ScintE.e4 = eta[3]>0. ? linearityCurve.applyLinCurve(3,t->ScintE.q4) / eta[3] : 0.;
    
    t->ScintE.nPE1 = eta[0] > 0. ? t->ScintE.e1 * eta[0] * alpha[0] : 0.;
    t->ScintE.nPE2 = eta[1] > 0. ? t->ScintE.e2 * eta[1] * alpha[1] : 0.;
    t->ScintE.nPE3 = eta[2] > 0. ? t->ScintE.e3 * eta[2] * alpha[2] : 0.;
    t->ScintE.nPE4 = eta[3] > 0. ? t->ScintE.e4 * eta[3] * alpha[3] : 0.;
    
    t->ScintE.de1 = t->ScintE.nPE1 > 0. ? t->ScintE.e1/sqrt(t->ScintE.nPE1) : 0.;
    t->ScintE.de2 = t->ScintE.nPE2 > 0. ? t->ScintE.e2/sqrt(t->ScintE.nPE2) : 0.;
    t->ScintE.de3 = t->ScintE.nPE3 > 0. ? t->ScintE.e3/sqrt(t->ScintE.nPE3) : 0.;
    t->ScintE.de4 = t->ScintE.nPE4 > 0. ? t->ScintE.e4/sqrt(t->ScintE.nPE4) : 0.;

    
    t->ScintW.e1 = eta[4]>0. ? linearityCurve.applyLinCurve(4,t->ScintW.q1) / eta[4] : 0.;
    t->ScintW.e2 = eta[5]>0. ? linearityCurve.applyLinCurve(5,t->ScintW.q2) / eta[5] : 0.;
    t->ScintW.e3 = eta[6]>0. ? linearityCurve.applyLinCurve(6,t->ScintW.q3) / eta[6] : 0.;
    t->ScintW.e4 = eta[7]>0. ? linearityCurve.applyLinCurve(7,t->ScintW.q4) / eta[7] : 0.;
    
    t->ScintW.nPE1 = eta[4] > 0. ? t->ScintW.e1 * eta[4] * alpha[4] : 0.;
    t->ScintW.nPE2 = eta[5] > 0. ? t->ScintW.e2 * eta[5] * alpha[5] : 0.;
    t->ScintW.nPE3 = eta[6] > 0. ? t->ScintW.e3 * eta[6] * alpha[6] : 0.;
    t->ScintW.nPE4 = eta[7] > 0. ? t->ScintW.e4 * eta[7] * alpha[7] : 0.;
    
    t->ScintW.de1 = t->ScintW.nPE1 > 0. ? t->ScintW.e1/sqrt(t->ScintW.nPE1) : 0.;
    t->ScintW.de2 = t->ScintW.nPE2 > 0. ? t->ScintW.e2/sqrt(t->ScintW.nPE2) : 0.;
    t->ScintW.de3 = t->ScintW.nPE3 > 0. ? t->ScintW.e3/sqrt(t->ScintW.nPE3) : 0.;
    t->ScintW.de4 = t->ScintW.nPE4 > 0. ? t->ScintW.e4/sqrt(t->ScintW.nPE4) : 0.;
    


    //Calculate the weighted energy on a side
    
    //EAST
    Double_t numer = (pmtQuality[0] ? t->ScintE.nPE1 : 0.) + (pmtQuality[1] ? t->ScintE.nPE2 : 0.) + (pmtQuality[2] ? t->ScintE.nPE3 : 0.) + (pmtQuality[3] ? t->ScintE.nPE4 : 0.);
    Double_t denom = 0.;
    for (int i=0; i<4; i++) denom += (pmtQuality[i] ? alpha[i] * eta[i] : 0.); 

    t->ScintE.energy = t->EvisE = (denom!=0. ? numer/denom : 0.);
    t->ScintE.denergy = (denom!=0. ? sqrt(t->ScintE.energy/denom) : 0.);

    //WEST
    numer = (pmtQuality[4] ? t->ScintW.nPE1 : 0.) + (pmtQuality[5] ? t->ScintW.nPE2 : 0.) + (pmtQuality[6] ? t->ScintW.nPE3 : 0.) + (pmtQuality[7] ? t->ScintW.nPE4 : 0.);
    denom = 0.;
    for (int i=4; i<8; i++) denom += (pmtQuality[i] ? alpha[i] * eta[i] : 0.); 

    t->ScintW.energy = t->EvisW = (denom!=0. ? numer/denom : 0.);
    t->ScintW.denergy = (denom!=0. ? sqrt(t->ScintW.energy/denom) : 0.);


    // Determine the reconstructed energy

    int typeIndex = t->Type==0 ? 0:(t->Type==1 ? 1:2); //for retrieving the parameters from EQ2Etrue
    double totalEvis=0.;
    
    if (t->Side==0) {
      totalEvis = t->Type==1 ? (t->EvisE+t->EvisW):t->EvisE;
      if (totalEvis>0.) {
	t->Erecon = EQ2Etrue[0][typeIndex][0]+EQ2Etrue[0][typeIndex][1]*totalEvis+EQ2Etrue[0][typeIndex][2]/(totalEvis+EQ2Etrue[0][typeIndex][3])+EQ2Etrue[0][typeIndex][4]/((totalEvis+EQ2Etrue[0][typeIndex][5])*(totalEvis+EQ2Etrue[0][typeIndex][5]));
      }
      else t->Erecon=-1.;
    }
    if (t->Side==1) {
      totalEvis = t->Type==1 ? (t->EvisE+t->EvisW):t->EvisW;
      if (totalEvis>0.) {
	t->Erecon = EQ2Etrue[1][typeIndex][0]+EQ2Etrue[1][typeIndex][1]*totalEvis+EQ2Etrue[1][typeIndex][2]/(totalEvis+EQ2Etrue[1][typeIndex][3])+EQ2Etrue[1][typeIndex][4]/((totalEvis+EQ2Etrue[1][typeIndex][5])*(totalEvis+EQ2Etrue[1][typeIndex][5]));
      }
      else t->Erecon=-1.;
    }
    
    t->fillOutputTree();
    
  }

  // Write output ntuple
  t->writeOutputFile();
  
  delete t; //Closes files

  return 0;
}

