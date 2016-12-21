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
#include "MBUtils.hh"
#include "positionMapHandler.hh"
#include "peaks.hh"
#include "calibrationTools.hh"

#include "replay_pass2.h"
#include "replay_pass3.h"


const bool useRCclasses = true;

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
  while (infile >> holdType >> params[side][type][0] >> params[side][type][1] >> params[side][type][2] >> params[side][type][3] 
	 >> params[side][type][4] >> params[side][type][5]) {
    std::cout << holdType << " " << params[side][type][0] << " " << params[side][type][1] << " " 
	      << params[side][type][2] << " " << params[side][type][3] << " " << params[side][type][4] 
	      << " " << params[side][type][5] << std::endl;
    type+=1;
    if (type==3) {type=0; side=1;}
  }
  return params;
}


vector <vector <double> > returnSourcePosition (Int_t runNumber, string src) {
  Char_t temp[500];
  sprintf(temp,"%s/source_list_%i.dat",getenv("SOURCE_LIST"),runNumber);
  ifstream file(temp);
  cout << src << endl;
  int num = 0;
  file >> num;
  cout << num << endl;
  int srcNum = 0;
  string src_check;
  for (int i=0; i<num;srcNum++,i++) {
    file >> src_check;
    cout << src_check << endl;
    if (src_check==src) break;   
  }
  cout << "The source Number is: " << srcNum << endl;
  if (srcNum==num) {
    cout << "Didn't find source in that run\n"; exit(0);
  }
  file.close();
  
  sprintf(temp,"%s/source_positions_%i.dat",getenv("SOURCE_POSITIONS"),runNumber);
  file.open(temp);
  
  vector < vector < double > > srcPos;
  srcPos.resize(2,vector <double> (3,0.));
  
  for (int i=0; i<srcNum+1; i++) {
    for (int j=0; j<2; j++) {
      for (int jj=0; jj<3; jj++) {
	file >> srcPos[j][jj];
      }
    }
  }
  return srcPos;
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
  //cout << "... Applying Calibration ..." << endl;

  unsigned int calibrationPeriod = getSrcRunPeriod(runNumber);
  // Determine linearity curve to use
  

  //Get the linearity curve
  LinearityCurve linearityCurve(calibrationPeriod, false);

  

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

  

  // Read source list
  int nSources = 0;
  string sourceName[3];
  char tempList[500];
  sprintf(tempList, "%s/source_list_%s.dat",getenv("SOURCE_LIST"), argv[1]);
  //cout << tempList << endl;
  ifstream fileList(tempList);
  if (fileList.is_open()) fileList >> nSources;
  cout << " ... Number of sources: " << nSources << endl;
  for (int n=0; n<nSources; n++) {
    fileList >> sourceName[n];
    cout << "  " << sourceName[n] << endl;
  }

  //Checking if one of the sources is Bi so we can add in that we will fit the second Bi peak as well
  bool useLowBiPeak=false;
  int BiPeakIndex = 0;
  for (int n=0; n<nSources; n++) {
    if (sourceName[n]=="Bi") {
      useLowBiPeak=true;
      BiPeakIndex=n;
      continue;
    }
  }
      

  //Adding in here that we are only looking for runs with Ce, Sn, In, or Bi in them
  bool correctSource = false;
  bool useSource[3] = {false,false,false};
  for (int n=0; n<nSources; n++) {
    if (sourceName[n]=="Ce" || sourceName[n]=="Sn" || sourceName[n]=="Bi" || sourceName[n]=="In") {
      correctSource = true;
      useSource[n]=true;
      continue; 
    }
  }
  
  if (!correctSource) {
    cout << "Source run with no Ce, Sn, or Bi\n";
    exit(0);
  }

  // Open output ntuple
  char tempOut[500];
  sprintf(tempOut, "%s/source_peaks_%s_Erecon.root",getenv("SOURCE_PEAKS"), argv[1]);
  TFile *fileOut = new TFile(tempOut,"RECREATE");


  // Read source positions
  double xEast[3], yEast[3], sigmaEast[3];
  double xWest[3], yWest[3], sigmaWest[3];
  char tempPos[500];
  sprintf(tempPos, "%s/source_positions_%s.dat",getenv("SOURCE_POSITIONS"), argv[1]);
  ifstream filePos(tempPos);
  cout << " ... Reading source positions" << endl;
  for (int n=0; n<nSources; n++) {
    filePos >> xEast[n] >> yEast[n] >> sigmaEast[n]
            >> xWest[n] >> yWest[n] >> sigmaWest[n];
    cout << xEast[n] << " " << yEast[n] << " " << sigmaEast[n] << " "
         << xWest[n] << " " << yWest[n] << " " << sigmaWest[n] << endl;
  }

  //Set up all the output histograms

  Int_t nBin = 600;
  TH1D *hisEvisTot[3][2];
  hisEvisTot[0][0] = new TH1D("hisEvis1E", "source1 East", nBin,-50.0,2400.0);
  hisEvisTot[0][1] = new TH1D("hisEvis1W", "source1 West", nBin,-50.0,2400.0);
  hisEvisTot[1][0] = new TH1D("hisEvis2E", "source2 East", nBin,-50.0,2400.0);
  hisEvisTot[1][1] = new TH1D("hisEvis2W", "source2 West", nBin,-50.0,2400.0);
  hisEvisTot[2][0] = new TH1D("hisEvis3E", "source3 East", nBin,-50.0,2400.0);
  hisEvisTot[2][1] = new TH1D("hisEvis3W", "source3 West", nBin,-50.0,2400.0);

  TH1D *hisEvis[3][8];
  hisEvis[0][0] = new TH1D("hisEvis1_E0", "", nBin,-50.0,2400.0);
  hisEvis[0][1] = new TH1D("hisEvis1_E1", "", nBin,-50.0,2400.0);
  hisEvis[0][2] = new TH1D("hisEvis1_E2", "", nBin,-50.0,2400.0);
  hisEvis[0][3] = new TH1D("hisEvis1_E3", "", nBin,-50.0,2400.0);
  hisEvis[0][4] = new TH1D("hisEvis1_W0", "", nBin,-50.0,2400.0);
  hisEvis[0][5] = new TH1D("hisEvis1_W1", "", nBin,-50.0,2400.0);
  hisEvis[0][6] = new TH1D("hisEvis1_W2", "", nBin,-50.0,2400.0);
  hisEvis[0][7] = new TH1D("hisEvis1_W3", "", nBin,-50.0,2400.0);

  hisEvis[1][0] = new TH1D("hisEvis2_E0", "", nBin,-50.0,2400.0);
  hisEvis[1][1] = new TH1D("hisEvis2_E1", "", nBin,-50.0,2400.0);
  hisEvis[1][2] = new TH1D("hisEvis2_E2", "", nBin,-50.0,2400.0);
  hisEvis[1][3] = new TH1D("hisEvis2_E3", "", nBin,-50.0,2400.0);
  hisEvis[1][4] = new TH1D("hisEvis2_W0", "", nBin,-50.0,2400.0);
  hisEvis[1][5] = new TH1D("hisEvis2_W1", "", nBin,-50.0,2400.0);
  hisEvis[1][6] = new TH1D("hisEvis2_W2", "", nBin,-50.0,2400.0);
  hisEvis[1][7] = new TH1D("hisEvis2_W3", "", nBin,-50.0,2400.0);

  hisEvis[2][0] = new TH1D("hisEvis3_E0", "", nBin,-50.0,2400.0);
  hisEvis[2][1] = new TH1D("hisEvis3_E1", "", nBin,-50.0,2400.0);
  hisEvis[2][2] = new TH1D("hisEvis3_E2", "", nBin,-50.0,2400.0);
  hisEvis[2][3] = new TH1D("hisEvis3_E3", "", nBin,-50.0,2400.0);
  hisEvis[2][4] = new TH1D("hisEvis3_W0", "", nBin,-50.0,2400.0);
  hisEvis[2][5] = new TH1D("hisEvis3_W1", "", nBin,-50.0,2400.0);
  hisEvis[2][6] = new TH1D("hisEvis3_W2", "", nBin,-50.0,2400.0);
  hisEvis[2][7] = new TH1D("hisEvis3_W3", "", nBin,-50.0,2400.0);
  

  TH1D *hisEreconTot[3][2];
  hisEreconTot[0][0] = new TH1D("hisErecon1E", "source1 East", nBin,0.0,2400.0);
  hisEreconTot[0][1] = new TH1D("hisErecon1W", "source1 West", nBin,0.0,2400.0);
  hisEreconTot[1][0] = new TH1D("hisErecon2E", "source2 East", nBin,0.0,2400.0);
  hisEreconTot[1][1] = new TH1D("hisErecon2W", "source2 West", nBin,0.0,2400.0);
  hisEreconTot[2][0] = new TH1D("hisErecon3E", "source3 East", nBin,0.0,2400.0);
  hisEreconTot[2][1] = new TH1D("hisErecon3W", "source3 West", nBin,0.0,2400.0);

  /*TH1D *hisErecon[3][8];
  his[0][0] = new TH1D("hisErecon1_E0", "", nBin,0.0,1200.0);
  his[0][1] = new TH1D("hisErecon1_E1", "", nBin,0.0,1200.0);
  his[0][2] = new TH1D("hisErecon1_E2", "", nBin,0.0,1200.0);
  his[0][3] = new TH1D("hisErecon1_E3", "", nBin,0.0,1200.0);
  his[0][4] = new TH1D("hisErecon1_W0", "", nBin,0.0,1200.0);
  his[0][5] = new TH1D("hisErecon1_W1", "", nBin,0.0,1200.0);
  his[0][6] = new TH1D("hisErecon1_W2", "", nBin,0.0,1200.0);
  his[0][7] = new TH1D("hisErecon1_W3", "", nBin,0.0,1200.0);

  his[1][0] = new TH1D("hisErecon2_E0", "", nBin,0.0,1200.0);
  his[1][1] = new TH1D("hisErecon2_E1", "", nBin,0.0,1200.0);
  his[1][2] = new TH1D("hisErecon2_E2", "", nBin,0.0,1200.0);
  his[1][3] = new TH1D("hisErecon2_E3", "", nBin,0.0,1200.0);
  his[1][4] = new TH1D("hisErecon2_W0", "", nBin,0.0,1200.0);
  his[1][5] = new TH1D("hisErecon2_W1", "", nBin,0.0,1200.0);
  his[1][6] = new TH1D("hisErecon2_W2", "", nBin,0.0,1200.0);
  his[1][7] = new TH1D("hisErecon2_W3", "", nBin,0.0,1200.0);

  his[2][0] = new TH1D("hisErecon3_E0", "", nBin,0.0,1200.0);
  his[2][1] = new TH1D("hisErecon3_E1", "", nBin,0.0,1200.0);
  his[2][2] = new TH1D("hisErecon3_E2", "", nBin,0.0,1200.0);
  his[2][3] = new TH1D("hisErecon3_E3", "", nBin,0.0,1200.0);
  his[2][4] = new TH1D("hisErecon3_W0", "", nBin,0.0,1200.0);
  his[2][5] = new TH1D("hisErecon3_W1", "", nBin,0.0,1200.0);
  his[2][6] = new TH1D("hisErecon3_W2", "", nBin,0.0,1200.0);
  his[2][7] = new TH1D("hisErecon3_W3", "", nBin,0.0,1200.0);*/
  
 
  // DataTree structure
  DataTree *t = new DataTree();

  // Open output ntuple
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

  vector < Double_t > eta;
  std::vector < std::vector <Int_t> > numDataPoints(3, std::vector <Int_t>(2,0)); //Holds the number of data points for each side for each source
  std::vector < std::vector <Double_t> > aveEta(3, std::vector <Double_t>(8,0.)); //Holds the average value of eta for the the data being read in for each source and each PMT
 
  
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


    ///////////////////////////////////////////////////////////////////////////
    // Now for filling the individual peak histograms using only type 0 events
    //////////////////////////////////////////////////////////////////////////
    if (t->Type != 0 || t->PID!=1) continue;

    double r2E = t->xE.center*t->xE.center + t->yE.center*t->yE.center;
    double r2W = t->xW.center*t->xW.center + t->yW.center*t->yW.center;
    double fiducialCut = 50.;

    if ( r2E>(fiducialCut*fiducialCut) || r2W>(fiducialCut*fiducialCut) ) continue; // get rid of events which are outside fiducial cut
    
    if (t->xE.nClipped>0 || t->yE.nClipped>0 || t->xW.nClipped>0 || t->yW.nClipped>0 ) continue; //Clipped events

    //  If the flag at the top of this file is set to true, also cut on the wirechamber
    //  event type according to C. Swanks classifications in ELOG 629 attachment 2
    if ( useRCclasses ) {
      if ( t->xeRC==0 || t->yeRC==0 || t->xwRC==0 || t->ywRC==0 ) continue;
    }

    if (useSource[0]) {
      
      if (t->Side == 0) {
	if ( (t->xE.center - xEast[0])*(t->xE.center - xEast[0]) +
	     (t->yE.center - yEast[0])*(t->yE.center - yEast[0]) <
	     (2.*sigmaEast[0])*(2.*sigmaEast[0]) ) {

	  numDataPoints[0][0]+=1;
	  aveEta[0][0] += eta[0];
	  aveEta[0][1] += eta[1];
	  aveEta[0][2] += eta[2];
	  aveEta[0][3] += eta[3];

	  if (t->Erecon>=0.) hisEreconTot[0][0]->Fill(t->Erecon);
	  hisEvisTot[0][0]->Fill(t->EvisE);

	  hisEvis[0][0]->Fill(t->ScintE.e1);
	  hisEvis[0][1]->Fill(t->ScintE.e2);
	  hisEvis[0][2]->Fill(t->ScintE.e3);
	  hisEvis[0][3]->Fill(t->ScintE.e4);
	  
	}
      }

      if (t->Side == 1) {
	if ( (t->xW.center - xWest[0])*(t->xW.center - xWest[0]) +
	     (t->yW.center - yWest[0])*(t->yW.center - yWest[0]) <
	     (2.*sigmaWest[0])*(2.*sigmaWest[0]) ) {

	  numDataPoints[0][1]+=1;
	  aveEta[0][4] += eta[4];
	  aveEta[0][5] += eta[5];
	  aveEta[0][6] += eta[6];
	  aveEta[0][7] += eta[7];
	  
	  if (t->Erecon>=0.) hisEreconTot[0][1]->Fill(t->Erecon);
	  hisEvisTot[0][1]->Fill(t->EvisW);

	  hisEvis[0][4]->Fill(t->ScintW.e1);
	  hisEvis[0][5]->Fill(t->ScintW.e2);
	  hisEvis[0][6]->Fill(t->ScintW.e3);
	  hisEvis[0][7]->Fill(t->ScintW.e4);

	}
      }
      
    }

    // Second source (x,y)
    if (nSources > 1 && useSource[1]) {


      if (t->Side == 0) {
	if ( (t->xE.center - xEast[1])*(t->xE.center - xEast[1]) +
	     (t->yE.center - yEast[1])*(t->yE.center - yEast[1]) <
	     (2.*sigmaEast[1])*(2.*sigmaEast[1]) ) {

	  numDataPoints[1][0]+=1;
	  aveEta[1][0] += eta[0];
	  aveEta[1][1] += eta[1];
	  aveEta[1][2] += eta[2];
	  aveEta[1][3] += eta[3];

	  if (t->Erecon>=0.) hisEreconTot[1][0]->Fill(t->Erecon);
	  hisEvisTot[1][0]->Fill(t->EvisE);

	  hisEvis[1][0]->Fill(t->ScintE.e1);
	  hisEvis[1][1]->Fill(t->ScintE.e2);
	  hisEvis[1][2]->Fill(t->ScintE.e3);
	  hisEvis[1][3]->Fill(t->ScintE.e4);

	}
      }
      if (t->Side == 1) {
	if ( (t->xW.center - xWest[1])*(t->xW.center - xWest[1]) +
	     (t->yW.center - yWest[1])*(t->yW.center - yWest[1]) <
	     (2.*sigmaWest[1])*(2.*sigmaWest[1]) ) {

	  numDataPoints[1][1]+=1;
	  aveEta[1][4] += eta[4];
	  aveEta[1][5] += eta[5];
	  aveEta[1][6] += eta[6];
	  aveEta[1][7] += eta[7];
	 
	  if (t->Erecon>=0.) hisEreconTot[1][1]->Fill(t->Erecon);
	  hisEvisTot[1][1]->Fill(t->EvisW);

	  hisEvis[1][4]->Fill(t->ScintW.e1);
	  hisEvis[1][5]->Fill(t->ScintW.e2);
	  hisEvis[1][6]->Fill(t->ScintW.e3);
	  hisEvis[1][7]->Fill(t->ScintW.e4);
 
	}
      }
    }

    // Third source (x,y)
    if (nSources > 2 && useSource[2]) {

      if (t->Side == 0) {
	if ( (t->xE.center - xEast[2])*(t->xE.center - xEast[2]) +
	     (t->yE.center - yEast[2])*(t->yE.center - yEast[2]) <
	     (2.*sigmaEast[2])*(2.*sigmaEast[2]) ) {

	  numDataPoints[2][0]+=1;
	  aveEta[2][0] += eta[0];
	  aveEta[2][1] += eta[1];
	  aveEta[2][2] += eta[2];
	  aveEta[2][3] += eta[3];
	  
	  if (t->Erecon>=0.) hisEreconTot[2][0]->Fill(t->Erecon);
	  hisEvisTot[2][0]->Fill(t->EvisE);

	  hisEvis[2][0]->Fill(t->ScintE.e1);
	  hisEvis[2][1]->Fill(t->ScintE.e2);
	  hisEvis[2][2]->Fill(t->ScintE.e3);
	  hisEvis[2][3]->Fill(t->ScintE.e4);
	}
      }
      if (t->Side == 1) {
	if ( (t->xW.center - xWest[2])*(t->xW.center - xWest[2]) +
	     (t->yW.center - yWest[2])*(t->yW.center - yWest[2]) <
	     (2.*sigmaWest[2])*(2.*sigmaWest[2]) ) {	  

	  numDataPoints[2][1]+=1;
	  aveEta[2][4] += eta[4];
	  aveEta[2][5] += eta[5];
	  aveEta[2][6] += eta[6];
	  aveEta[2][7] += eta[7];
	   
	  if (t->Erecon>=0.) hisEreconTot[2][1]->Fill(t->Erecon);
	  hisEvisTot[2][1]->Fill(t->EvisW);

	  hisEvis[2][4]->Fill(t->ScintW.e1);
	  hisEvis[2][5]->Fill(t->ScintW.e2);
	  hisEvis[2][6]->Fill(t->ScintW.e3);
	  hisEvis[2][7]->Fill(t->ScintW.e4);

	}
      }

    }
    
  }

  // Write output ntuple
  t->writeOutputFile();
  
  delete t; //Closes files

  //Finish average eta calculation and write to file
  {
    std::string etaFile = getenv("SOURCE_POSITIONS")+std::string("meanEta_")+itos(runNumber)+".dat";
    ofstream meanEtaVals(etaFile.c_str());
    for (int src=0; src<3; src++) {
      for (int pmt=0; pmt<8; pmt++) {
	
	aveEta[src][pmt] = useSource[src] ? (pmt<4 ? aveEta[src][pmt]/(double)numDataPoints[src][0] : aveEta[src][pmt]/(double)numDataPoints[src][1]) : 0.;
	
      }
    }

    for (int n=0; n<3; n++) {
      meanEtaVals << runNumber << " "
		  << sourceName[n] << " "
		  << aveEta[n][0] << " "
		  << aveEta[n][1] << " "
		  << aveEta[n][2] << " "
		  << aveEta[n][3] << " "
		  << aveEta[n][4] << " "
		  << aveEta[n][5] << " "
		  << aveEta[n][6] << " "
		  << aveEta[n][7] << endl;
    }
    meanEtaVals.close();
  }


  
  // Find maximum bin
  double maxBin[3][8]={0.};
  double maxCounts[3][8]={0.};
  for (int n=0; n<nSources; n++) {
    for (int j=0; j<8; j++) {
      
      maxBin[n][j] = hisEvis[n][j]->GetMaximumBin();
      maxCounts[n][j] = hisEvis[n][j]->GetBinContent(maxBin[n][j]);

    }
  }

  // Define histogram fit ranges
  double xLow[3][8]={0.}, xHigh[3][8]={0.};
  for (int n=0; n<nSources; n++) {
    
    for (int j=0; j<8; j++) {
      for (int i=maxBin[n][j]; i<nBin; i++) {
        if (hisEvis[n][j]->GetBinContent(i+1) < 0.25*maxCounts[n][j]) {
          xHigh[n][j] = hisEvis[n][j]->GetBinCenter(i+1);
	  //Check to make sure the value isn't too close to the maximum bin...
          if ((i-maxBin[n][j])<5) xHigh[n][j] = hisEvis[n][j]->GetXaxis()->GetBinCenter(maxBin[n][j])+250.;
	  break;
        }
	if (i==(nBin-1)) xHigh[n][j] = hisEvis[n][j]->GetBinCenter(i);
      }
      if (sourceName[n]!="Bi") {
	for (int i=maxBin[n][j]; i>0; i--) {
	  if (hisEvis[n][j]->GetBinContent(i-1) < 0.25*maxCounts[n][j]) {
	    xLow[n][j] = hisEvis[n][j]->GetBinCenter(i-1);
	    //if ((maxBin[n][j]-i)<5) xLow[n][j] = binCenterMax[n][j]-200.;
	    break;
	  }
	}
      }
      else {
	for (int i=maxBin[n][j]*0.5; i>0; i--) {
	  if (hisEvis[n][j]->GetBinContent(i-1) < 0.2*0.5*maxCounts[n][j]) {
	    xLow[n][j] = hisEvis[n][j]->GetBinCenter(i-1);
	    //if ((maxBin[n][j]-i)<5) xLow[n][j] = binCenterMax[n][j]-200.;
	    break;
	  }
	}
      }
    }
    
  }

  double fitMean[3][8]={0.};
  double fitMeanError[3][8]={0.};
  double lowBiFitMean[8]={0.};
  double lowBiFitMeanError[8]={0.};
  double fitSigma[3][8]={0.};
  double lowBiFitSigma[8]={0.};

  for (int n=0; n<nSources; n++) {

    if (useSource[n]) {

      for (int j=0; j<8; j++) {

	if (sourceName[n]!="Bi") {

	  //Double_t rangeLow = 5.;
	  //Double_t rangeHigh = 4096.;
	  SinglePeakHist sing(hisEvis[n][j], xLow[n][j], xHigh[n][j]);

	  if (sing.isGoodFit()) {
	    fitMean[n][j] = sing.ReturnMean();
	    fitMeanError[n][j] = sing.ReturnMeanError();
	    fitSigma[n][j] = sing.ReturnSigma();
	  }

	  else  {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " peak in PMT " << j << ". Trying one more time......" << endl;
	    sing.SetRangeMin(xLow[n][j]);
	    sing.SetRangeMax(xHigh[n][j]);
	    sing.FitHist(maxBin[n][j], 40., hisEvis[n][j]->GetBinContent(maxBin[n][j]));

	    if (sing.isGoodFit()) { 
	      fitMean[n][j] = sing.ReturnMean();
	      fitMeanError[n][j] = sing.ReturnMeanError();
	      fitSigma[n][j] = sing.ReturnSigma();
	    }

	    else cout << "RUN " << runNumber << " CAN'T CONVERGE ON " << sourceName[n] << " PEAK IN PMT " << j  << endl;
	  }
	}

	else {

	  //Double_t rangeLow = 5.;
	  //Double_t rangeHigh = 4096.;
	  DoublePeakHist doub(hisEvis[n][j], xLow[n][j], xHigh[n][j]);

	  if (doub.isGoodFit()) {
	    //std::cout << "IT WAS A GOOD FIT\n";
	    fitMean[n][j] = doub.ReturnMean1();
	    fitMeanError[n][j] = doub.ReturnMean1Error();
	    lowBiFitMean[j] = doub.ReturnMean2();
	    lowBiFitMeanError[j] = doub.ReturnMean2Error();
	    fitSigma[n][j] = doub.ReturnSigma1();
	    lowBiFitSigma[j] = doub.ReturnSigma2();
	  }

	  else  {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " peak in PMT " << j << ". Trying one more time......" << endl;
	    doub.SetRangeMin(xLow[n][j]);
	    doub.SetRangeMax(xHigh[n][j]);
	    doub.FitHist(maxBin[n][j], 40., hisEvis[n][j]->GetBinContent(maxBin[n][j]), 0.5*maxBin[n][j], 40., 0.5*hisEvis[n][j]->GetBinContent(maxBin[n][j]));

	    if (doub.isGoodFit()) {
	      fitMean[n][j] = doub.ReturnMean1();
	      fitMeanError[n][j] = doub.ReturnMean1Error();
	      lowBiFitMean[j] = doub.ReturnMean2();
	      lowBiFitMeanError[j] = doub.ReturnMean2Error();
	      fitSigma[n][j] = doub.ReturnSigma1();
	      lowBiFitSigma[j] = doub.ReturnSigma2();
	    }

	    else cout << "RUN " << runNumber << " CAN'T CONVERGE ON " << sourceName[n] << " PEAK IN PMT " << j  << endl;
	    
	  }
	}

      }
    }
  }

  // Write results to file
  char tempResults[500];
  sprintf(tempResults, "%s/source_peaks_%s_Evis.dat",getenv("SOURCE_PEAKS"), argv[1]);
  ofstream outResultsMean(tempResults);
  sprintf(tempResults, "%s/source_peaks_errors_%s_Evis.dat",getenv("SOURCE_PEAKS"), argv[1]);
  ofstream outResultsMeanError(tempResults);
  sprintf(tempResults, "%s/source_widths_%s_Evis.dat",getenv("SOURCE_PEAKS"), argv[1]);
  ofstream outResultsSigma(tempResults);

  //sprintf(tempResults, "%s/source_peaks_%s_Evis.dat","likelihoodCheck", argv[1]);
  //ofstream outResultsMean(tempResults);
  //sprintf(tempResults, "%s/source_peaks_errors_%s_Evis.dat","likelihoodCheck", argv[1]);
  //ofstream outResultsMeanError(tempResults);
  //sprintf(tempResults, "%s/source_widths_%s_Evis.dat","likelihoodCheck", argv[1]);
  //ofstream outResultsSigma(tempResults);

  for (int n=0; n<nSources; n++) {
    if (useSource[n]) {
      string src = sourceName[n];
      if (sourceName[n]=="Bi") src=sourceName[n]+"1";
      outResultsMean << runNumber << " "
		     << src << " "
		     << fitMean[n][0] << " "
		     << fitMean[n][1] << " "
		     << fitMean[n][2] << " "
		     << fitMean[n][3] << " "
		     << fitMean[n][4] << " "
		     << fitMean[n][5] << " "
		     << fitMean[n][6] << " "
		     << fitMean[n][7] << endl;
      outResultsMeanError << runNumber << " "
		     << src << " "
		     << fitMeanError[n][0] << " "
		     << fitMeanError[n][1] << " "
		     << fitMeanError[n][2] << " "
		     << fitMeanError[n][3] << " "
		     << fitMeanError[n][4] << " "
		     << fitMeanError[n][5] << " "
		     << fitMeanError[n][6] << " "
		     << fitMeanError[n][7] << endl;
      outResultsSigma << runNumber << " "
		     << src << " "
		     << fitSigma[n][0] << " "
		     << fitSigma[n][1] << " "
		     << fitSigma[n][2] << " "
		     << fitSigma[n][3] << " "
		     << fitSigma[n][4] << " "
		     << fitSigma[n][5] << " "
		     << fitSigma[n][6] << " "
		     << fitSigma[n][7] << endl;
    }
  }

  if (useLowBiPeak) {
    outResultsMean << runNumber << " "
		   << "Bi2" << " "
		   << lowBiFitMean[0] << " "
		   << lowBiFitMean[1] << " "
		   << lowBiFitMean[2] << " "
		   << lowBiFitMean[3] << " "
		   << lowBiFitMean[4] << " "
		   << lowBiFitMean[5] << " "
		   << lowBiFitMean[6] << " "
		   << lowBiFitMean[7] << endl;
    outResultsMeanError << runNumber << " "
		   << "Bi2" << " "
		   << lowBiFitMeanError[0] << " "
		   << lowBiFitMeanError[1] << " "
		   << lowBiFitMeanError[2] << " "
		   << lowBiFitMeanError[3] << " "
		   << lowBiFitMeanError[4] << " "
		   << lowBiFitMeanError[5] << " "
		   << lowBiFitMeanError[6] << " "
		   << lowBiFitMeanError[7] << endl;
    outResultsSigma << runNumber << " "
		   << "Bi2" << " "
		   << lowBiFitSigma[0] << " "
		   << lowBiFitSigma[1] << " "
		   << lowBiFitSigma[2] << " "
		   << lowBiFitSigma[3] << " "
		   << lowBiFitSigma[4] << " "
		   << lowBiFitSigma[5] << " "
		   << lowBiFitSigma[6] << " "
		   << lowBiFitSigma[7] << endl;
  }
  outResultsMean.close();
  outResultsMeanError.close();
  outResultsSigma.close();


  //Now for the total Evis peaks and widths

  // Find maximum bin
  double maxBinEvisTot[3][2]={0.};
  double maxCountsEvisTot[3][2]={0.};
  for (int n=0; n<nSources; n++) {
    for (int j=0; j<2; j++) {
      maxBinEvisTot[n][j] = hisEvisTot[n][j]->GetMaximumBin();
      maxCountsEvisTot[n][j] = hisEvisTot[n][j]->GetBinContent(maxBinEvisTot[n][j]);
    }
  }


  // Define histogram fit ranges
  double xLowEvisTot[3][2]={0.}, xHighEvisTot[3][2]={0.};
  for (int n=0; n<nSources; n++) { 
    for (int j=0; j<2; j++) {
      for (int i=maxBinEvisTot[n][j]; i<nBin; i++) {
	if (hisEvisTot[n][j]->GetBinContent(i+1) < 0.33*maxCountsEvisTot[n][j]) {
	  xHighEvisTot[n][j] = hisEvisTot[n][j]->GetBinCenter(i+1);
	  //Check to make sure the value isn't too close to the maximum bin...
	  if ((i-maxBinEvisTot[n][j])<5) xHighEvisTot[n][j] = hisEvisTot[n][j]->GetXaxis()->GetBinCenter(maxBinEvisTot[n][j])+250.;
	  break;
	}
	if (i==(nBin-1)) xHighEvisTot[n][j] = hisEvisTot[n][j]->GetBinCenter(i);
      }
      if (sourceName[n].substr(0,2)!="Bi") {
	for (int i=maxBinEvisTot[n][j]; i>0; i--) {
	  if (hisEvisTot[n][j]->GetBinContent(i-1) < 0.33*maxCountsEvisTot[n][j]) {
	    xLowEvisTot[n][j] = hisEvisTot[n][j]->GetBinCenter(i-1);
	    //if ((maxBinEvisTot[n][j]-i)<5) xLowEvisTot[n][j] = binCenterMax[n][j]-200.;
	    break;
	  }
	}
      }
      else {
	for (int i=maxBinEvisTot[n][j]*0.5; i>0; i--) {
	  if (hisEvisTot[n][j]->GetBinContent(i-1) < 0.33*0.5*maxCountsEvisTot[n][j]) {
	    xLowEvisTot[n][j] = hisEvisTot[n][j]->GetBinCenter(i-1);
	    //if ((maxBinEvisTot[n][j]-i)<5) xLowEvisTot[n][j] = binCenterMax[n][j]-200.;
	    break;
	  }
	}
      }
    }
  }

  double fitMeanEvisTot[3][2]={0.};
  double lowBiFitMeanEvisTot[2]={0.};
  double fitSigmaEvisTot[3][2]={0.};
  double lowBiFitSigmaEvisTot[2]={0.};

  for (int n=0; n<nSources; n++) {

    if (useSource[n]) {
      
      for (int j=0; j<2; j++) {

	if (sourceName[n].substr(0,2)!="Bi") {
	
	  
	  //Double_t rangeLowEvisTot = 5.;
	  //Double_t rangeHighEvisTot = 4096.;
	  SinglePeakHist sing(hisEvisTot[n][j], xLowEvisTot[n][j], xHighEvisTot[n][j]);
	  
	  if (sing.isGoodFit()) {
	    fitMeanEvisTot[n][j] = sing.ReturnMean();
	    fitSigmaEvisTot[n][j] = sing.ReturnSigma();
	  }
	  
	  else  {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " EvisTot peak.... Trying one more time......" << endl;
	    sing.SetRangeMin(xLowEvisTot[n][j]);
	    sing.SetRangeMax(xHighEvisTot[n][j]);
	    sing.FitHist(maxBinEvisTot[n][j], 40., hisEvisTot[n][j]->GetBinContent(maxBinEvisTot[n][j]));
	    
	    if (sing.isGoodFit()) { 
	      fitMeanEvisTot[n][j] = sing.ReturnMean();
	      fitSigmaEvisTot[n][j] = sing.ReturnSigma();
	    }
	    
	    else cout << "RUN " << runNumber << " CAN'T CONVERGE ON " << sourceName[n] << " EvisTot PEAK "<<  endl;
	  }
	}
	
	else {
	  
	  
	  //DoublePeakHist doub(hisEvisTot[n][j], xLowEvisTot[n][j], xHighEvisTot[n][j]);
	  SinglePeakHist singBi1(hisEvisTot[n][j], 800., xHighEvisTot[n][j]);
	  SinglePeakHist singBi2(hisEvisTot[n][j], 340., 530.);
	  
	  if (singBi1.isGoodFit()) {	    
	    fitMeanEvisTot[n][j] = singBi1.ReturnMean();	 
	    fitSigmaEvisTot[n][j] = singBi1.ReturnSigma();
	  }
	  else  {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " EvisTot peak.... Trying one more time......" << endl;
	    singBi1.SetRangeMin(775.);
	    singBi1.SetRangeMax(xHighEvisTot[n][j]);
	    singBi1.FitHist(910., 40., hisEvisTot[n][j]->GetBinContent(maxBinEvisTot[n][j]));
	    
	    if (singBi1.isGoodFit()) {
	      fitMeanEvisTot[n][j] = singBi1.ReturnMean();
	      fitSigmaEvisTot[n][j] = singBi1.ReturnSigma();
	    }
	    
	    else cout << "RUN " << runNumber << " CAN'T CONVERGE ON " << sourceName[n] << " EvisTot PEAK " << endl;
	    
	  }
	  
	  if (singBi2.isGoodFit()) {	    
	    lowBiFitMeanEvisTot[j] = singBi2.ReturnMean();	 
	    lowBiFitSigmaEvisTot[j] = singBi2.ReturnSigma();
	  }
	  else  {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " EvisTot peak.... Trying one more time......" << endl;
	    singBi2.SetRangeMin(340.);
	    singBi2.SetRangeMax(530.);
	    singBi2.FitHist(430., 40., 0.4*hisEvisTot[n][j]->GetBinContent(maxBinEvisTot[n][j]));
	    
	    if (singBi2.isGoodFit()) {
	      lowBiFitMeanEvisTot[j] = singBi2.ReturnMean();
	      lowBiFitSigmaEvisTot[j] = singBi2.ReturnSigma();
	    }
	    
	    else cout << "RUN " << runNumber << " CAN'T CONVERGE ON " << sourceName[n] << " EvisTot PEAK " << endl;
	  }
	}
      }
      
    }
  }
  
  
  sprintf(tempResults, "%s/source_peaks_%s_EvisTot.dat",getenv("SOURCE_PEAKS"), argv[1]);
  ofstream outResultsEvisTot(tempResults);
  
  for (int n=0; n<nSources; n++) {
    if (useSource[n]) {
      string src = sourceName[n];
      if (sourceName[n]=="Bi") src=sourceName[n]+"1";
     
      outResultsEvisTot << runNumber << " "
		       << src << " "
		       << fitMeanEvisTot[n][0] << " "
		       << fitSigmaEvisTot[n][0] << " "
		       << fitMeanEvisTot[n][1] << " " 
		       << fitSigmaEvisTot[n][1] << std::endl;
      
    }
  }
  
  if (useLowBiPeak) {
    outResultsEvisTot << runNumber << " "
		     << "Bi2" << " "
		     << lowBiFitMeanEvisTot[0] << " "
		     << lowBiFitSigmaEvisTot[0] << " "
		     << lowBiFitMeanEvisTot[1] << " " 
		     << lowBiFitSigmaEvisTot[1] << std::endl;
  }
  outResultsEvisTot.close();




  //Now for the Erecon peaks and widths

  // Find maximum bin
  double maxBinEreconTot[3][2]={0.};
  double maxCountsEreconTot[3][2]={0.};
  for (int n=0; n<nSources; n++) {
    for (int j=0; j<2; j++) {
      maxBinEreconTot[n][j] = hisEreconTot[n][j]->GetMaximumBin();
      maxCountsEreconTot[n][j] = hisEreconTot[n][j]->GetBinContent(maxBinEreconTot[n][j]);
    }
  }


  // Define histogram fit ranges
  double xLowEreconTot[3][2]={0.}, xHighEreconTot[3][2]={0.};
  for (int n=0; n<nSources; n++) { 
    for (int j=0; j<2; j++) {
      for (int i=maxBinEreconTot[n][j]; i<nBin; i++) {
	if (hisEreconTot[n][j]->GetBinContent(i+1) < 0.33*maxCountsEreconTot[n][j]) {
	  xHighEreconTot[n][j] = hisEreconTot[n][j]->GetBinCenter(i+1);
	  //Check to make sure the value isn't too close to the maximum bin...
	  if ((i-maxBinEreconTot[n][j])<5) xHighEreconTot[n][j] = hisEreconTot[n][j]->GetXaxis()->GetBinCenter(maxBinEreconTot[n][j])+250.;
	  break;
	}
	if (i==(nBin-1)) xHighEreconTot[n][j] = hisEreconTot[n][j]->GetBinCenter(i);
      }
      if (sourceName[n].substr(0,2)!="Bi") {
	for (int i=maxBinEreconTot[n][j]; i>0; i--) {
	  if (hisEreconTot[n][j]->GetBinContent(i-1) < 0.33*maxCountsEreconTot[n][j]) {
	    xLowEreconTot[n][j] = hisEreconTot[n][j]->GetBinCenter(i-1);
	    //if ((maxBinEreconTot[n][j]-i)<5) xLowEreconTot[n][j] = binCenterMax[n][j]-200.;
	    break;
	  }
	}
      }
      else {
	for (int i=maxBinEreconTot[n][j]*0.5; i>0; i--) {
	  if (hisEreconTot[n][j]->GetBinContent(i-1) < 0.33*0.5*maxCountsEreconTot[n][j]) {
	    xLowEreconTot[n][j] = hisEreconTot[n][j]->GetBinCenter(i-1);
	    //if ((maxBinEreconTot[n][j]-i)<5) xLowEreconTot[n][j] = binCenterMax[n][j]-200.;
	    break;
	  }
	}
      }
    }
  }

  double fitMeanEreconTot[3][2]={0.};
  double lowBiFitMeanEreconTot[2]={0.};
  double fitSigmaEreconTot[3][2]={0.};
  double lowBiFitSigmaEreconTot[2]={0.};

  for (int n=0; n<nSources; n++) {

    if (useSource[n]) {
      
      for (int j=0; j<2; j++) {

	if (sourceName[n].substr(0,2)!="Bi") {
	
	  
	  //Double_t rangeLowEreconTot = 5.;
	  //Double_t rangeHighEreconTot = 4096.;
	  SinglePeakHist sing(hisEreconTot[n][j], xLowEreconTot[n][j], xHighEreconTot[n][j]);
	  
	  if (sing.isGoodFit()) {
	    fitMeanEreconTot[n][j] = sing.ReturnMean();
	    fitSigmaEreconTot[n][j] = sing.ReturnSigma();
	  }
	  
	  else  {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " EreconTot peak.... Trying one more time......" << endl;
	    sing.SetRangeMin(xLowEreconTot[n][j]);
	    sing.SetRangeMax(xHighEreconTot[n][j]);
	    sing.FitHist(maxBinEreconTot[n][j], 40., hisEreconTot[n][j]->GetBinContent(maxBinEreconTot[n][j]));
	    
	    if (sing.isGoodFit()) { 
	      fitMeanEreconTot[n][j] = sing.ReturnMean();
	      fitSigmaEreconTot[n][j] = sing.ReturnSigma();
	    }
	    
	    else cout << "RUN " << runNumber << " CAN'T CONVERGE ON " << sourceName[n] << " EreconTot PEAK "<<  endl;
	  }
	}
	
	else {
	  
	  
	  //DoublePeakHist doub(hisEreconTot[n][j], xLowEreconTot[n][j], xHighEreconTot[n][j]);
	  SinglePeakHist singBi1(hisEreconTot[n][j], 850., xHighEreconTot[n][j]);
	  SinglePeakHist singBi2(hisEreconTot[n][j], 390., 580.);
	  
	  if (singBi1.isGoodFit()) {	    
	    fitMeanEreconTot[n][j] = singBi1.ReturnMean();	 
	    fitSigmaEreconTot[n][j] = singBi1.ReturnSigma();
	  }
	  else  {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " EreconTot peak.... Trying one more time......" << endl;
	    singBi1.SetRangeMin(850.);
	    singBi1.SetRangeMax(990.);
	    singBi1.FitHist(960., 40., hisEreconTot[n][j]->GetBinContent(maxBinEreconTot[n][j]));
	    
	    if (singBi1.isGoodFit()) {
	      fitMeanEreconTot[n][j] = singBi1.ReturnMean();
	      fitSigmaEreconTot[n][j] = singBi1.ReturnSigma();
	    }
	    
	    else cout << "RUN " << runNumber << " CAN'T CONVERGE ON " << sourceName[n] << " EreconTot PEAK " << endl;
	    
	  }
	  
	  if (singBi2.isGoodFit()) {	    
	    lowBiFitMeanEreconTot[j] = singBi2.ReturnMean();	 
	    lowBiFitSigmaEreconTot[j] = singBi2.ReturnSigma();
	  }
	  else  {
	    cout << "Run " << runNumber << " can't converge on " << sourceName[n] << " EreconTot peak.... Trying one more time......" << endl;
	    singBi2.SetRangeMin(390.);
	    singBi2.SetRangeMax(580.);
	    singBi2.FitHist(480., 40., 0.4*hisEreconTot[n][j]->GetBinContent(maxBinEreconTot[n][j]));
	    
	    if (singBi2.isGoodFit()) {
	      lowBiFitMeanEreconTot[j] = singBi2.ReturnMean();
	      lowBiFitSigmaEreconTot[j] = singBi2.ReturnSigma();
	    }
	    
	    else cout << "RUN " << runNumber << " CAN'T CONVERGE ON " << sourceName[n] << " EreconTot PEAK " << endl;
	  }
	}
      }
      
    }
  }
  
  
  sprintf(tempResults, "%s/source_peaks_%s_EreconTot.dat",getenv("SOURCE_PEAKS"), argv[1]);
  ofstream outResultsEreconTot(tempResults);
  
  for (int n=0; n<nSources; n++) {
    if (useSource[n]) {
      string src = sourceName[n];
      if (sourceName[n]=="Bi") src=sourceName[n]+"1";
       
      outResultsEreconTot << runNumber << " "
		       << src << " "
		       << fitMeanEreconTot[n][0] << " "
		       << fitSigmaEreconTot[n][0] << " "
		       << fitMeanEreconTot[n][1] << " " 
		       << fitSigmaEreconTot[n][1] << std::endl;
      
    }
  }
  
  if (useLowBiPeak) {
    outResultsEreconTot << runNumber << " "
		     << "Bi2" << " "
		     << lowBiFitMeanEreconTot[0] << " "
		     << lowBiFitSigmaEreconTot[0] << " "
		     << lowBiFitMeanEreconTot[1] << " " 
		     << lowBiFitSigmaEreconTot[1] << std::endl;
  }
  outResultsEreconTot.close();
  

  // NOW I CALCULATE THE ADJUSTED ADC VALUE TO USE FOR CALIBRATION PURPOSES
  //std::vector < std::vector < Double_t > > ADCbar(nSources, std::vector <Double_t> (8,0.));
  //for (int n=0; n<nSources; n++) {
  // for (int p=0; p<8; p++) {
  //  ADCbar[n][p] = linearityCurve.applyInverseLinCurve(p,fitMean[n][p]*aveEta[n][p<4?0:1]);
  //  std::cout << "Source " << n << " PMT " << p << " ADC :" << aveEta[n][p] << std::endl;
  //}
  //}
 
  sprintf(tempResults, "%s/source_peaks_%s_ADC.dat",getenv("SOURCE_PEAKS"), argv[1]);
  ofstream outResultsADC(tempResults);
  sprintf(tempResults, "%s/source_peaks_errors_%s_ADC.dat",getenv("SOURCE_PEAKS"), argv[1]);
  ofstream outResultsADCError(tempResults);

  std::vector < std::vector < Double_t > > calParams = linearityCurve.returnLinCurveParams();

  for (int n=0; n<nSources; n++) {
    if (useSource[n]) {
      std::vector < std::vector <Double_t> > pos = returnSourcePosition(runNumber, sourceName[n]); //Holds the average value of eta for the the data being read in for each source and each PMT
      std::vector < Double_t > eta0 = posmap.getInterpolatedEta(pos[0][0], pos[0][1], pos[1][0], pos[1][1]);

      for (Int_t ii=0; ii<8; ii++) {
	//if (aveEta[n][ii] > 2.5 || aveEta[n][ii] < 0.25) aveEta[n][ii] = eta0[ii];
	aveEta[n][ii] = eta0[ii];
	std::cout << "Source " << sourceName[n] << " PMT " << ii << " aveEta :" << aveEta[n][ii] << std::endl;
      }

      string src = sourceName[n];
      if (sourceName[n]=="Bi") src=sourceName[n]+"1";
  
      outResultsADC << runNumber << " "
		    << src << " "
		    << linearityCurve.applyInverseLinCurve(0,fitMean[n][0]*aveEta[n][0]) << " "
		    << linearityCurve.applyInverseLinCurve(1,fitMean[n][1]*aveEta[n][1]) << " "
		    << linearityCurve.applyInverseLinCurve(2,fitMean[n][2]*aveEta[n][2]) << " "
		    << linearityCurve.applyInverseLinCurve(3,fitMean[n][3]*aveEta[n][3]) << " "
		    << linearityCurve.applyInverseLinCurve(4,fitMean[n][4]*aveEta[n][4]) << " "
		    << linearityCurve.applyInverseLinCurve(5,fitMean[n][5]*aveEta[n][5]) << " "
		    << linearityCurve.applyInverseLinCurve(6,fitMean[n][6]*aveEta[n][6]) << " "
		    << linearityCurve.applyInverseLinCurve(7,fitMean[n][7]*aveEta[n][7]) << endl;
      outResultsADCError << runNumber << " "
			 << src << " "
			 << fitMeanError[n][0]/calParams[0][1] << " "
			 << fitMeanError[n][1]/calParams[1][1] << " "
			 << fitMeanError[n][2]/calParams[2][1] << " "
			 << fitMeanError[n][3]/calParams[3][1] << " "
			 << fitMeanError[n][4]/calParams[4][1] << " "
			 << fitMeanError[n][5]/calParams[5][1] << " "
			 << fitMeanError[n][6]/calParams[6][1] << " "
			 << fitMeanError[n][7]/calParams[7][1] << " " << endl;
    }
  }

  if (useLowBiPeak) {

    std::vector < std::vector <Double_t> > pos = returnSourcePosition(runNumber, "Bi"); //Holds the average value of eta for the the data being read in for each source and each PMT
    std::vector < Double_t > eta0 = posmap.getInterpolatedEta(pos[0][0], pos[0][1], pos[1][0], pos[1][1]);
    
    for (Int_t ii=0; ii<8; ii++) {
      aveEta[BiPeakIndex][ii] = eta0[ii];
      std::cout << "Source " << "Bi2" << " PMT " << ii << " aveEta :" << aveEta[BiPeakIndex][ii] << std::endl;
    }

    outResultsADC << runNumber << " "
		  << "Bi2" << " "
		  << linearityCurve.applyInverseLinCurve(0,lowBiFitMean[0]*aveEta[BiPeakIndex][0]) << " "
		  << linearityCurve.applyInverseLinCurve(1,lowBiFitMean[1]*aveEta[BiPeakIndex][1]) << " "
		  << linearityCurve.applyInverseLinCurve(2,lowBiFitMean[2]*aveEta[BiPeakIndex][2]) << " "
		  << linearityCurve.applyInverseLinCurve(3,lowBiFitMean[3]*aveEta[BiPeakIndex][3]) << " "
		  << linearityCurve.applyInverseLinCurve(4,lowBiFitMean[4]*aveEta[BiPeakIndex][4]) << " "
		  << linearityCurve.applyInverseLinCurve(5,lowBiFitMean[5]*aveEta[BiPeakIndex][5]) << " "
		  << linearityCurve.applyInverseLinCurve(6,lowBiFitMean[6]*aveEta[BiPeakIndex][6]) << " "
		  << linearityCurve.applyInverseLinCurve(7,lowBiFitMean[7]*aveEta[BiPeakIndex][7]) << std::endl;
		
    outResultsADCError << runNumber << " "
		       << "Bi2" << " "
		       << lowBiFitMeanError[0]/calParams[0][1] << " "
		       << lowBiFitMeanError[1]/calParams[1][1] << " "
		       << lowBiFitMeanError[2]/calParams[2][1] << " "
		       << lowBiFitMeanError[3]/calParams[3][1] << " "
		       << lowBiFitMeanError[4]/calParams[4][1] << " "
		       << lowBiFitMeanError[5]/calParams[5][1] << " "
		       << lowBiFitMeanError[6]/calParams[6][1] << " "
		       << lowBiFitMeanError[7]/calParams[7][1] << " " << endl;
  }
  outResultsADC.close();

  // Write output ntuple
  fileOut->Write();
  fileOut->Close();

  return 0;
}




