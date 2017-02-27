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
#include "MWPCPositionResponse.hh"



#include "replay_pass2.h"
#include "replay_pass3.h"

bool OnlyReplayBadFiles = false;

std::vector < std::vector <Double_t> > loadPMTpedestals(Int_t runNumber) {

  Char_t temp[500];
  std::vector < std::vector < Double_t > > peds (8,std::vector <Double_t> (2,0.));
  sprintf(temp,"%s/PMT_pedestals_%i.dat",getenv("PEDESTALS"),runNumber);
  ifstream infile;
  infile.open(temp);

  Int_t i = 0;
  Int_t run;

  while (infile >> run >> peds[i][0] >> peds[i][1]) { std::cout << "Pedestal " << i << ": " << peds[i][0] << " " << peds[i][1] << std::endl; i++; }
  return peds;

};


//Get the conversion from EQ2Etrue
std::vector < std::vector < std::vector <double> > > getEQ2EtrueParams(int runNumber) {
  ifstream infile;
  std::string basePath = getenv("ANALYSIS_CODE"); 
  if (runNumber<16000) basePath+=std::string("/simulation_comparison/EQ2EtrueConversion/2011-2012_EQ2EtrueFitParams.dat");
  else if (runNumber<20000) basePath+=std::string("/simulation_comparison/EQ2EtrueConversion/2011-2012_EQ2EtrueFitParams.dat");
  else if (runNumber<21628 && runNumber>21087) basePath+=std::string("/simulation_comparison/EQ2EtrueConversion/2012-2013_isobutane_EQ2EtrueFitParams.dat");
  else if (runNumber<24000) basePath+=std::string("/simulation_comparison/EQ2EtrueConversion/2012-2013_EQ2EtrueFitParams.dat");
  else {
    std::cout << "Bad runNumber passed to getEQ2EtrueParams\n";
    exit(0);
  }
  infile.open(basePath.c_str());
  std::vector < std::vector < std::vector < double > > > params;
  params.resize(2,std::vector < std::vector < double > > (3, std::vector < double > (6,0.)));

  char holdType[10];
  int side=0, type=0;
  while (infile >> holdType >> params[side][type][0] >> params[side][type][1] >> params[side][type][2] >> params[side][type][3] >> params[side][type][4] >> params[side][type][5]) { 
    std::cout << holdType << " " << params[side][type][0] << " " << params[side][type][1] << " " << params[side][type][2] << " " << params[side][type][3] << " " << params[side][type][4] << " " << params[side][type][5] << std::endl;
    type+=1;
    if (type==3) {type=0; side=1;}
  }
  return params;
}

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
  int nParams = 3; //takes a quadratic 

  // Run number integer
  cout << "Run " << runNumber << " ..." << endl;

  char tempOut[500];
  sprintf(tempOut, "%s/replay_pass3_%s.root",getenv("REPLAY_PASS3"), argv[1]);
  //sprintf(tempOut, "replay_pass3_%s.root", argv[1]);

  //Check if the file is good already and quit if it is so that we can 
  // only replay the files that are bad...

  if ( OnlyReplayBadFiles ) {
     
    if ( checkIfReplayFileIsGood(std::string(tempOut)) == 1 ) return 1;
  
    else {
      std::ofstream badRuns("badRuns.txt", std::fstream::app);
      badRuns << argv[1] << "\n";
      badRuns.close();
    }
  }

  // Reading in pedestals file to get cathode pedestals
  // Read pedestals file
  char tempFilePed[500];
  int iRun;
  sprintf(tempFilePed, "%s/pedestals_%s.dat", getenv("PEDESTALS"), argv[1]);
  std::cout << "... Reading: " << tempFilePed << std::endl;

  std::ifstream filePed(tempFilePed);
  for (int i=0; i<8; i++) {
    filePed >> iRun >> pedQadc[i];   
  }
  for (int i=0; i<32; i++) {
    filePed >> iRun >> pedPdc2[i];
  }
  for (int i=0; i<32; i++) {
    filePed >> iRun >> pedPadc[i];
  }

  
  filePed >> iRun >> pedPdc30;
  filePed >> iRun >> pedPdc34;
  


  cout << "... Applying Calibration ..." << endl;

  unsigned int calibrationPeriod = getSrcRunPeriod(runNumber);

  //Get the linearity curve
  LinearityCurve linearityCurve(calibrationPeriod,false);

  //Load the simulated relationship between EQ and Etrue
  vector < vector < vector < double > > > EQ2Etrue = getEQ2EtrueParams(runNumber);
 
  //Read in PMT quality file
  std::vector <Int_t> pmtQuality = getPMTQuality(runNumber);

  //Get values for nPE/keV...
  std::vector <Double_t> alpha = GetAlphaValues(calibrationPeriod);
  
  //Reading position map...
  UInt_t XePeriod = getXeRunPeriod(runNumber); // Get the proper Xe run period for the Trigger functions
  //GetPositionMap(XePeriod);
  PositionMap posmap(5.0,50.);
  posmap.readPositionMap(XePeriod);

  // DataTree structure
  DataTree *t = new DataTree();

  // Open output ntuple
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

    if ( i%10000==0 ) std::cout << i << std::endl;

    // Only process event if it's an electron
    
    if ( t->PID==1 ) {

      std::vector <double> posex(3,0.);
      std::vector <double> poswx(3,0.);
      std::vector <double> posey(3,0.);
      std::vector <double> poswy(3,0.);

      MWPCCathodeHandler cathResp(t->Cathodes_Ex,t->Cathodes_Ey,t->Cathodes_Wx,t->Cathodes_Wy,&pedPdc2[16],&pedPdc2[0],&pedPadc[16],&pedPadc[0]);
      
      cathResp.findAllPositions();

      posex = cathResp.getPosEX();
      posey = cathResp.getPosEY();
      poswx = cathResp.getPosWX();
      poswy = cathResp.getPosWY();

      t->xE.center = posex[0] * positionProjection;
      t->yE.center = posey[0] * positionProjection;
      t->xW.center = poswx[0] * positionProjection;
      t->yW.center = poswy[0] * positionProjection;

      t->xE.width = posex[1] * positionProjection;
      t->yE.width = posey[1] * positionProjection;
      t->xW.width = poswx[1] * positionProjection;
      t->yW.width = poswy[1] * positionProjection;
      
      t->xE.height = posex[2];
      t->yE.height = posey[2];
      t->xW.height = poswx[2];
      t->yW.height = poswy[2];

      t->xE.mult = cathResp.getMultEX();
      t->yE.mult = cathResp.getMultEY();
      t->xW.mult = cathResp.getMultWX();
      t->yW.mult = cathResp.getMultWY();

      t->xE.nClipped = cathResp.getnClippedEX();
      t->yE.nClipped = cathResp.getnClippedEY();
      t->xW.nClipped = cathResp.getnClippedWX();
      t->yW.nClipped = cathResp.getnClippedWY();

      t->xE.maxWire = cathResp.getMaxWireEX();
      t->yE.maxWire = cathResp.getMaxWireEY();
      t->xW.maxWire = cathResp.getMaxWireWX();
      t->yW.maxWire = cathResp.getMaxWireWY();

      t->xE.maxValue = t->Cathodes_Ex[t->xE.maxWire];
      t->yE.maxValue = t->Cathodes_Ey[t->yE.maxWire];
      t->xW.maxValue = t->Cathodes_Wx[t->xW.maxWire];
      t->yW.maxValue = t->Cathodes_Wy[t->yW.maxWire];

      t->xE.rawCenter = cathResp.getWirePosEX(t->xE.maxWire);
      t->yE.rawCenter = cathResp.getWirePosEY(t->yE.maxWire);
      t->xW.rawCenter = cathResp.getWirePosWX(t->xW.maxWire);
      t->yW.rawCenter = cathResp.getWirePosWY(t->yW.maxWire);

      // write out pedestal subtracted cathode values
      
      for ( int ii = 0; ii<16; ++ii ) {
	t->Cathodes_Ex[ii] = t->Cathodes_Ex[ii] - pedPdc2[ii+16]; 
	t->Cathodes_Ey[ii] = t->Cathodes_Ey[ii] - pedPdc2[ii];
	t->Cathodes_Wx[ii] = t->Cathodes_Wx[ii] - pedPadc[ii+16];
	t->Cathodes_Wy[ii] = t->Cathodes_Wy[ii] - pedPadc[ii];
      }

      /////// Now do the energy reconstruction
      
      eta = posmap.getInterpolatedEta(t->xE.center, t->yE.center, t->xW.center, t->yW.center);
      
      t->ScintE.e1 = linearityCurve.applyLinCurve(0,t->ScintE.q1);
      t->ScintE.e2 = linearityCurve.applyLinCurve(1,t->ScintE.q2);
      t->ScintE.e3 = linearityCurve.applyLinCurve(2,t->ScintE.q3);
      t->ScintE.e4 = linearityCurve.applyLinCurve(3,t->ScintE.q4);
      
      
      t->ScintE.e1 = ( eta[0]>0. && t->ScintE.e1>0. ) ? t->ScintE.e1 / eta[0] : 0.;
      t->ScintE.e2 = ( eta[1]>0. && t->ScintE.e2>0. ) ? t->ScintE.e2 / eta[1] : 0.;
      t->ScintE.e3 = ( eta[2]>0. && t->ScintE.e3>0. ) ? t->ScintE.e3 / eta[2] : 0.;
      t->ScintE.e4 = ( eta[3]>0. && t->ScintE.e4>0. ) ? t->ScintE.e4 / eta[3] : 0.;
      
      t->ScintE.nPE1 = eta[0] > 0. ? t->ScintE.e1 * eta[0] * alpha[0] : 0.;
      t->ScintE.nPE2 = eta[1] > 0. ? t->ScintE.e2 * eta[1] * alpha[1] : 0.;
      t->ScintE.nPE3 = eta[2] > 0. ? t->ScintE.e3 * eta[2] * alpha[2] : 0.;
      t->ScintE.nPE4 = eta[3] > 0. ? t->ScintE.e4 * eta[3] * alpha[3] : 0.;
      
      t->ScintE.de1 = t->ScintE.nPE1 > 0. ? t->ScintE.e1/sqrt(t->ScintE.nPE1) : 0.;
      t->ScintE.de2 = t->ScintE.nPE2 > 0. ? t->ScintE.e2/sqrt(t->ScintE.nPE2) : 0.;
      t->ScintE.de3 = t->ScintE.nPE3 > 0. ? t->ScintE.e3/sqrt(t->ScintE.nPE3) : 0.;
      t->ScintE.de4 = t->ScintE.nPE4 > 0. ? t->ScintE.e4/sqrt(t->ScintE.nPE4) : 0.;
      
      
      t->ScintW.e1 = linearityCurve.applyLinCurve(4,t->ScintW.q1);
      t->ScintW.e2 = linearityCurve.applyLinCurve(5,t->ScintW.q2);
      t->ScintW.e3 = linearityCurve.applyLinCurve(6,t->ScintW.q3);
      t->ScintW.e4 = linearityCurve.applyLinCurve(7,t->ScintW.q4);
      
      t->ScintW.e1 = ( eta[4]>0. && t->ScintW.e1>0. ) ? t->ScintW.e1 / eta[4] : 0.;
      t->ScintW.e2 = ( eta[5]>0. && t->ScintW.e2>0. ) ? t->ScintW.e2 / eta[5] : 0.;
      t->ScintW.e3 = ( eta[6]>0. && t->ScintW.e3>0. ) ? t->ScintW.e3 / eta[6] : 0.;
      t->ScintW.e4 = ( eta[7]>0. && t->ScintW.e4>0. ) ? t->ScintW.e4 / eta[7] : 0.;
      
      t->ScintW.nPE1 = eta[4] > 0. ? t->ScintW.e1 * eta[4] * alpha[4] : 0.;
      t->ScintW.nPE2 = eta[5] > 0. ? t->ScintW.e2 * eta[5] * alpha[5] : 0.;
      t->ScintW.nPE3 = eta[6] > 0. ? t->ScintW.e3 * eta[6] * alpha[6] : 0.;
      t->ScintW.nPE4 = eta[7] > 0. ? t->ScintW.e4 * eta[7] * alpha[7] : 0.;
      
      t->ScintW.de1 = t->ScintW.nPE1 > 0. ? t->ScintW.e1/sqrt(t->ScintW.nPE1) : 0.;
      t->ScintW.de2 = t->ScintW.nPE2 > 0. ? t->ScintW.e2/sqrt(t->ScintW.nPE2) : 0.;
      t->ScintW.de3 = t->ScintW.nPE3 > 0. ? t->ScintW.e3/sqrt(t->ScintW.nPE3) : 0.;
      t->ScintW.de4 = t->ScintW.nPE4 > 0. ? t->ScintW.e4/sqrt(t->ScintW.nPE4) : 0.;
      
      //std::cout << "Made it here" << std::endl;
      
      
      //Calculate the weighted energy on a side
      
      //EAST
      Double_t numer = ( (pmtQuality[0] && t->ScintE.nPE1>0. ? t->ScintE.nPE1 : 0.) +
			 (pmtQuality[1] && t->ScintE.nPE1>0. ? t->ScintE.nPE2 : 0.) + 
			 (pmtQuality[2] && t->ScintE.nPE1>0. ? t->ScintE.nPE3 : 0.) + 
			 (pmtQuality[3] && t->ScintE.nPE1>0. ? t->ScintE.nPE4 : 0.) );
      
      Double_t denom  = ( (pmtQuality[0] && t->ScintE.nPE1>0. ? alpha[0] * eta[0] : 0.) +
			  (pmtQuality[1] && t->ScintE.nPE2>0. ? alpha[1] * eta[1] : 0.) +
			  (pmtQuality[2] && t->ScintE.nPE3>0. ? alpha[2] * eta[2] : 0.) + 
			  (pmtQuality[3] && t->ScintE.nPE4>0. ? alpha[3] * eta[3] : 0.) ); 
      
      t->ScintE.energy = t->EvisE = (denom!=0. ? numer/denom : 0.);
      t->ScintE.denergy = (denom!=0. ? sqrt(t->ScintE.energy/denom) : 0.);
      
      //WEST
      numer = denom = 0.;
      
      numer = ( (pmtQuality[4] && t->ScintW.nPE1>0. ? t->ScintW.nPE1 : 0.) +
		(pmtQuality[5] && t->ScintW.nPE1>0. ? t->ScintW.nPE2 : 0.) + 
		(pmtQuality[6] && t->ScintW.nPE1>0. ? t->ScintW.nPE3 : 0.) + 
		(pmtQuality[7] && t->ScintW.nPE1>0. ? t->ScintW.nPE4 : 0.) );
      
      denom  = ( (pmtQuality[4] && t->ScintW.nPE1>0. ? alpha[4] * eta[4] : 0.) +
		 (pmtQuality[5] && t->ScintW.nPE2>0. ? alpha[5] * eta[5] : 0.) +
		 (pmtQuality[6] && t->ScintW.nPE3>0. ? alpha[6] * eta[6] : 0.) + 
		 (pmtQuality[7] && t->ScintW.nPE4>0. ? alpha[7] * eta[7] : 0.) ); 
      
      
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
    }
    
    t->fillOutputTree();
    
  }
  // Write output ntuple
  t->writeOutputFile();
  
  delete t; //Closes files

  return 0;
}

