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

#include "replay_pass2.h"
#include "replay_pass3.h"

using namespace std;

int main(int argc, char *argv[])
{
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(12);

  // Position bins
  int nPMT = 8;
  int nPosBinsX = 43;
  int nPosBinsY = 43;  
  double xBinWidth = 2.5;
  double yBinWidth = 2.5;
  double xBinLower[nPosBinsX];
  double xBinUpper[nPosBinsX];
  double xBinCenter[nPosBinsX];
  double yBinLower[nPosBinsY];
  double yBinUpper[nPosBinsY];
  double yBinCenter[nPosBinsY];
  int intXBinCenter[nPosBinsX];
  int intYBinCenter[nPosBinsY];

  for (int k=0; k<nPosBinsX; k++) {
    xBinLower[k]     = -(double)nPosBinsX*xBinWidth/2. + ((double) k)*xBinWidth;
    xBinUpper[k]     = -(double)nPosBinsX*xBinWidth/2. + ((double) k)*xBinWidth + xBinWidth;
    xBinCenter[k]    = (xBinLower[k] + xBinUpper[k])/2.;
    intXBinCenter[k] = (int) xBinCenter[k];
    //cout << xBinLower[k] << " " << intXBinCenter[k] << " " << xBinUpper[k] << endl;
  }

  for (int k=0; k<nPosBinsY; k++) {
    yBinLower[k]     = -(double)nPosBinsY*yBinWidth/2. + ((double) k)*yBinWidth;
    yBinUpper[k]     = -(double)nPosBinsY*yBinWidth/2. + ((double) k)*yBinWidth + yBinWidth;
    yBinCenter[k]    = (yBinLower[k] + yBinUpper[k])/2.;
    intYBinCenter[k] = (int) yBinCenter[k];
    //cout << yBinLower[k] << " " << intYBinCenter[k] << " " << yBinUpper[k] << endl;
  }

  

  // Run number integer
  cout << "Run " << argv[1] << " ..." << endl;
  cout << "... Applying Xe position map ..." << endl;
  istringstream ss(argv[1]);
  int runNumber;
  ss >> runNumber;

  // Determine position map to use
  char tempFileXePositionMap[500];

  unsigned int XePeriod = getXeRunPeriod(runNumber);
  sprintf(tempFileXePositionMap, "../position_map/position_map_%i_RC_123.dat",XePeriod);
  
  cout << "... Reading: " << tempFileXePositionMap << endl;

  // Read position map
  double x, y;
  double positionMap[nPMT][nPosBinsX][nPosBinsY];
  ifstream fileXePositionMap(tempFileXePositionMap);
  for (int i=0; i<nPosBinsX; i++) {
    for (int j=0; j<nPosBinsY; j++) {
      fileXePositionMap >> x >> y
                        >> positionMap[0][i][j]
			>> positionMap[1][i][j]
			>> positionMap[2][i][j]
			>> positionMap[3][i][j]
			>> positionMap[4][i][j]
			>> positionMap[5][i][j]
			>> positionMap[6][i][j]
			>> positionMap[7][i][j];
    }
  }

  // Calculate (x,y) correction factor: eta = 1/positionMap
  double eta[nPMT][nPosBinsX][nPosBinsY];
  for (int i=0; i<nPosBinsX; i++) {
    for (int j=0; j<nPosBinsY; j++) {
      eta[0][i][j] = 1. / positionMap[0][i][j];
      eta[1][i][j] = 1. / positionMap[1][i][j];
      eta[2][i][j] = 1. / positionMap[2][i][j];
      eta[3][i][j] = 1. / positionMap[3][i][j];
      eta[4][i][j] = 1. / positionMap[4][i][j];
      eta[5][i][j] = 1. / positionMap[5][i][j];
      eta[6][i][j] = 1. / positionMap[6][i][j];
      eta[7][i][j] = 1. / positionMap[7][i][j];
    }
  }

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

  // Loop over events
  for (int i=0; i<nEvents; i++) {
    t->getEvent(i);

    // Determine (x,y) bin
    int intEastBinX = -1;
    int intEastBinY = -1;
    int intWestBinX = -1;
    int intWestBinY = -1;

    for (int m=0; m<nPosBinsX; m++) {
      if ( (t->xE.center >= xBinLower[m]) && (t->xE.center < xBinUpper[m]) ) intEastBinX = m;
      if ( (t->xW.center >= xBinLower[m]) && (t->xW.center < xBinUpper[m]) ) intWestBinX = m;
    }

    for (int m=0; m<nPosBinsY; m++) {
      if ( (t->yE.center >= yBinLower[m]) && (t->yE.center < yBinUpper[m]) ) intEastBinY = m;
      if ( (t->yW.center >= yBinLower[m]) && (t->yW.center < yBinUpper[m]) ) intWestBinY = m;
    }

    // Apply (x,y) correction factor
    if (intEastBinX > -1 && intEastBinY > -1) {
      t->ScintE.q1 = t->ScintE.q1 * eta[0][intEastBinX][intEastBinY];
      t->ScintE.q2 = t->ScintE.q2 * eta[1][intEastBinX][intEastBinY];
      t->ScintE.q3 = t->ScintE.q3 * eta[2][intEastBinX][intEastBinY];
      t->ScintE.q4 = t->ScintE.q4 * eta[3][intEastBinX][intEastBinY];
    }
    
    if (intWestBinX > -1 && intWestBinY > -1) {
      t->ScintW.q1 = t->ScintW.q1 * eta[4][intWestBinX][intWestBinY];
      t->ScintW.q2 = t->ScintW.q2 * eta[5][intWestBinX][intWestBinY];
      t->ScintW.q3 = t->ScintW.q3 * eta[6][intWestBinX][intWestBinY];
      t->ScintW.q4 = t->ScintW.q4 * eta[7][intWestBinX][intWestBinY];
    }

    
    t->fillOutputTree();
  }

  // Write output ntuple
  t->writeOutputFile();
  
  delete t; //Closes files

  return 0;
}

//////////// OLD INPUT OUTPUT FILE METHOD /////////////////
/*
TFile *fileOut = new TFile(tempOut,"RECREATE");
  TTree *Tout = new TTree("pass3", "pass3");

  // Variables
  Tout->Branch("pmt0_pass3", &pmt_pass3[0], "pmt0_pass3/D");
  Tout->Branch("pmt1_pass3", &pmt_pass3[1], "pmt1_pass3/D");
  Tout->Branch("pmt2_pass3", &pmt_pass3[2], "pmt2_pass3/D");
  Tout->Branch("pmt3_pass3", &pmt_pass3[3], "pmt3_pass3/D");
  Tout->Branch("pmt4_pass3", &pmt_pass3[4], "pmt4_pass3/D");
  Tout->Branch("pmt5_pass3", &pmt_pass3[5], "pmt5_pass3/D");
  Tout->Branch("pmt6_pass3", &pmt_pass3[6], "pmt6_pass3/D");
  Tout->Branch("pmt7_pass3", &pmt_pass3[7], "pmt7_pass3/D");

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

  Tout->Branch("xE_pass3", &xE_pass3, "xE_pass3/D");
  Tout->Branch("yE_pass3", &yE_pass3, "yE_pass3/D");
  Tout->Branch("xW_pass3", &xW_pass3, "xW_pass3/D");
  Tout->Branch("yW_pass3", &yW_pass3, "yW_pass3/D");
  
  int xeRC,yeRC,xwRC,ywRC;	
  Tout->Branch("xeRC", &xeRC, "xeRC/I"); //x east response class. 
  Tout->Branch("yeRC", &yeRC, "yeRC/I"); //y east response class... 
  Tout->Branch("xwRC", &xwRC, "xwRC/I");
  Tout->Branch("ywRC", &ywRC, "ywRC/I");

  Tout->Branch("PID_pass3",  &PID_pass3,  "PID_pass3/I");
  Tout->Branch("type_pass3", &type_pass3, "type_pass3/I");
  Tout->Branch("side_pass3", &side_pass3, "side_pass3/I");
  Tout->Branch("posError_pass3", &posError_pass3, "posError_pass3/I");

  // Open input ntuple
  char tempIn[500];
  sprintf(tempIn, "%s/replay_pass2_%s.root", getenv("REPLAY_PASS2"),argv[1]);

  TFile *fileIn = new TFile(tempIn, "READ");
  TTree *Tin = (TTree*)(fileIn->Get("pass2"));

  // Variables
  Tin->SetBranchAddress("pmt0_pass2", &pmt_pass2[0]);
  Tin->SetBranchAddress("pmt1_pass2", &pmt_pass2[1]);
  Tin->SetBranchAddress("pmt2_pass2", &pmt_pass2[2]);
  Tin->SetBranchAddress("pmt3_pass2", &pmt_pass2[3]);
  Tin->SetBranchAddress("pmt4_pass2", &pmt_pass2[4]);
  Tin->SetBranchAddress("pmt5_pass2", &pmt_pass2[5]);
  Tin->SetBranchAddress("pmt6_pass2", &pmt_pass2[6]);
  Tin->SetBranchAddress("pmt7_pass2", &pmt_pass2[7]);

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

  Tin->SetBranchAddress("xE_pass2", &xE_pass2);
  Tin->SetBranchAddress("yE_pass2", &yE_pass2);
  Tin->SetBranchAddress("xW_pass2", &xW_pass2);
  Tin->SetBranchAddress("yW_pass2", &yW_pass2);

  //Adding in wirechamber class variable
  Tin->SetBranchAddress("xeRC",&xeRC);
  Tin->SetBranchAddress("yeRC",&yeRC);
  Tin->SetBranchAddress("xwRC",&xwRC);
  Tin->SetBranchAddress("ywRC",&ywRC);

  Tin->SetBranchAddress("PID_pass2",  &PID_pass2);
  Tin->SetBranchAddress("type_pass2", &type_pass2);
  Tin->SetBranchAddress("side_pass2", &side_pass2);
  Tin->SetBranchAddress("posError_pass2", &posError_pass2);
*/
